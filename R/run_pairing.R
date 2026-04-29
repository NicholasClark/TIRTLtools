#' Find TCRalpha/beta pairs from individual well read counts
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This runs the MAD-HYPE and T-SHELL algorithms to find TCRalpha-beta pairs
#' originating from the same clone.
#'
#'
#' @param folder_path the path of the folder with well-level data
#' @param folder_out the path of the folder to write results to. The function will create the folder if it does not exist.
#' @param prefix a prefix for the output file names
#' @param well_filter_thres wells are removed if they have fewer unique clones than: wellfilter_thres*(Avg. # of unique clones per well).
#' The default value is 0.5
#' @param min_reads minimum number of reads a chain must have in a well to be considered observed
#' (note: actual minimum is min_reads+1. Default value is 0, i.e. chain must have >= 1 read in a well)
#' @param min_wells minimum number of wells a chain must be observed in to be paired.
#' @param well_pos the position of the well ID (e.g. "B5") in the file names. For example, files
#' named "<lab>_<project>_<well_id>_TCRalpha.tsv" would use well_pos=3. (default is 3)
#' @param wellset1 a vector of wells to use for the pairing
#' @param compute whether or not to run the pairing algorithms after tabulating and writing pseudobulk data (default TRUE)
#' @param backend the computing backend to use. The function looks for a GPU and automatically chooses an appropriate backend by default.
#' @param pval_thres_tshell the adjusted p-value threshold for T-SHELL significance (default 1e-10)
#' @param wij_thres_tshell the threshold for the number of wells containing both chains for T-SHELL significance (default >2 wells)
#' @param verbose whether to print out messages (default TRUE)
#' @param write_extra_files whether to write un-necessary intermediate files (default FALSE)
#' @param filter_before_top3 whether to filter by loss fraction before extracting top 3 correlation values for T-SHELL (default FALSE)
#' @param fork whether to "fork" the python process for basilisk (default is NULL, which automatically chooses an appropriate option)
#' @param shared whether to use a "shared" python process for basilisk (default is NULL, which automatically chooses an appropriate option)
#' @param chunk_size batch size for calculations in pairing scripts
#' @param exclude_nonfunctional whether to exclude non-functional chains before pairing (default is FALSE)
#' @param select_best_madhype whether to use a secondary algorithm on the pairs from the MAD-HYPE algorithm to select
#' the best pairs for each clone (default is FALSE)
#' @param select_best_tshell whether to use a secondary algorithm on the pairs from the T-SHELL algorithm to select
#' the best pairs for each clone (default is FALSE)
#'
#' @return
#' A data frame with the TCR-alpha/TCR-beta pairs.
#'
#' The function also writes three files to "folder_out":
#'  - A data frame ("<prefix>_pseudobulk_TRA.tsv") of pseudobulk counts for TCRalpha chains
#'  - A data frame ("<prefix>_pseudobulk_TRB.tsv") of pseudobulk counts for TCRbeta chains
#'  - A data frame ("<prefix>_TIRTLoutput.tsv") of TCR-alpha/TCR-beta pairs.
#'
#'  These files can be loaded using the \code{\link{load_tirtlseq}()} function.
#'
#'  If write_extra_files is TRUE, the function also writes sparse matrices of per-well read
#'  counts (well x chain) for TCR-alpha and beta to "<prefix>_alpha_mat.rds" and "<prefix>_beta_mat.rds".
#'  Metadata for the chains in these matrices are written to "<prefix>_alpha_meta.parquet" and
#'  "<prefix>_beta_meta.parquet" and metadata for the wells is written to "<prefix>_well_meta.parquet".
#'
#'  These files can be loaded using the \code{\link{load_well_counts_binary}()} function.
#'
#' @family pairing
#' @export
#'

run_pairing = function(
    folder_path,
    folder_out,
    prefix,
    well_filter_thres=0.5,
    min_reads=0,
    min_wells=2,
    well_pos=3,
    wellset1=get_well_subset(1:16,1:24),
    compute=TRUE,
    backend = c("auto", "cpu", "cupy", "mlx"),
    pval_thres_tshell=1e-10,
    wij_thres_tshell=2,
    verbose = TRUE,
    write_extra_files = FALSE,
    filter_before_top3 = FALSE,
    fork = NULL,
    shared = NULL,
    chunk_size = 500,
    exclude_nonfunctional = FALSE,
    select_best_madhype = FALSE,
    select_best_tshell = FALSE
){
  tictoc::tic()
  #ensure_python_env()
  py_require( packages = get_py_deps() )

  backend = backend[1]

  .create_folder(folder_out)

  df_alpha = tibble(file_alpha = list.files(path = folder_path,full.names = F, pattern = ".*_TRA\\.tsv")) %>%
    mutate(file_alpha_no_dot = gsub(".","_",file_alpha,fixed=T),
           well = strsplit(file_alpha_no_dot, "_") %>% sapply(., function(x) x[[well_pos]]),
           well_row = substr(well, 0, 1),
           well_column = substr(well,2,3) %>% as.integer(),
           in_wellset = well %in% wellset1) %>%
    filter(in_wellset) %>%
    arrange(well_row, well_column) %>%
    select(-well_row, -well_column, -in_wellset)
  df_beta = tibble(file_beta = list.files(path = folder_path,full.names = F, pattern = ".*_TRB\\.tsv")) %>%
    mutate(file_beta_no_dot = gsub(".","_",file_beta,fixed=T),
           well = strsplit(file_beta_no_dot, "_") %>% sapply(., function(x) x[[well_pos]]),
           well_row = substr(well, 0, 1),
           well_column = substr(well,2,3) %>% as.integer(),
           in_wellset = well %in% wellset1)  %>%
    filter(in_wellset) %>%
    arrange(well_row, well_column) %>%
    select(-well_row, -well_column, -in_wellset)
  missing_alpha_wells = wellset1[!wellset1 %in% df_alpha$well]
  missing_beta_wells = wellset1[!wellset1 %in% df_beta$well]
  n_files_alpha = nrow(df_alpha)
  n_files_beta = nrow(df_beta)
  if(length(missing_alpha_wells) > 0) {
    msg = paste(length(missing_alpha_wells), "wells missing .tsv file (TCRalpha):", paste0(missing_alpha_wells, collapse = ", "))
    if(verbose) message(msg)
  }
  if(length(missing_beta_wells) > 0) {
    msg = paste(length(missing_beta_wells), "wells missing .tsv file (TCRbeta):", paste0(missing_beta_wells, collapse = ", "))
    if(verbose) message(msg)
  }
  df_meta = full_join(df_alpha, df_beta, by = "well")
  wells_missing = wellset1[!wellset1 %in% df_meta$well]
  df_meta = df_meta %>%
    bind_rows(tibble(well = wells_missing)) %>%
    mutate(well_row = substr(well, 0, 1), well_column = substr(well,2,3) %>% as.integer()) %>%
    arrange(well_row, well_column)
  # missing_wells_final = wellset1[!wellset1 %in% df_meta$well]
  # if(length(missing_wells_final) > 0) {
  #   msg = paste(length(missing_wells_final), "wells excluded for missing .tsv files:", paste0(missing_wells_final, collapse = ", "))
  #   if(verbose) message(msg)
  # }

  # files = list.files(path = folder_path,full.names = F)
  # files_no_dot = gsub(".","_",files,fixed=T)
  # file_wells = strsplit(files_no_dot, "_") %>% sapply(., function(x) x[[well_pos]])
  # files_bool = file_wells %in% wellset1
  # msg = paste("Reading", sum(files_bool), "files from", folder_path)
  # if(verbose) message(msg)
  # files_load = file.path(folder_path, files[files_bool])

  if(verbose) message(paste("Loading clone files (TCRalpha) for", n_files_alpha, "wells...", collapse=" "))
  n_files_loaded_alpha = 0
  mlista = lapply(1:nrow(df_meta),function(i){
    ff = df_meta$file_alpha[i]
    if(!is.na(ff)) {
      tmp = fread(file.path(folder_path, ff))
      n_files_loaded_alpha <<- n_files_loaded_alpha + 1
      if(verbose) if(n_files_loaded_alpha %% 25 == 0 | n_files_loaded_alpha == n_files_alpha) message(paste(n_files_loaded_alpha, "of", n_files_alpha, "files loaded"))
    } else {
      tmp = tibble()
    }
    return(tmp)
  }) %>% set_names(df_meta$well)
  if(verbose) message(paste("Loading clone files (TCRbeta) for", n_files_beta, "wells...", collapse=" "))
  n_files_loaded_beta = 0
  mlistb = lapply(1:nrow(df_meta),function(i){
    ff = df_meta$file_beta[i]
    if(!is.na(ff)) {
      tmp = fread(file.path(folder_path, ff))
      n_files_loaded_beta <<- n_files_loaded_beta + 1
      if(verbose) if(n_files_loaded_beta %% 25 == 0 | n_files_loaded_beta == n_files_beta) message(paste(n_files_loaded_beta, "of", n_files_beta, "files loaded"))
    } else {
      tmp = tibble()
    }
    return(tmp)
  }) %>% set_names(df_meta$well)

  # mlist<-lapply(1:length(files_load),function(i){
  #   if(verbose) if(i %% 25 == 0 | i == length(files_load)) message(paste(i, "of", length(files_load), "files loaded"))
  #   fread(files_load[i])
  # })
  # names(mlist) = files_no_dot[files_bool]

  # mlista<-geta(mlist)
  # mlistb<-getb(mlist)
  # if(verbose) {
  #   n_filesA = length(mlista)
  #   n_filesB = length(mlistb)
  #   msgA = paste(n_filesA, "TCRalpha well files loaded")
  #   msgB = paste(n_filesB, "TCRbeta well files loaded")
  #   message(msgA)
  #   message(msgB)
  # }
  tmpa = bind_rows(mlista)
  tmpb = bind_rows(mlistb)
  n_counts_char_alpha = sum(tmpa$readCount) %>% scales::scientific()
  n_counts_char_beta = sum(tmpb$readCount) %>% scales::scientific()
  msg1 = paste("Total number of TCRalpha reads:", n_counts_char_alpha)
  msg2 = paste("Total number of TCRbeta reads:", n_counts_char_beta)
  if(verbose) message(msg1)
  if(verbose) message(msg2)
  n_unique_alpha = length(unique(tmpa$targetSequences))
  n_unique_beta = length(unique(tmpb$targetSequences))
  msg1 = paste("Total number of unique TCRalpha chains:", n_unique_alpha)
  msg2 = paste("Total number of unique TCRbeta chains:", n_unique_beta)
  if(verbose) message(msg1)
  if(verbose) message(msg2)

  alpha_nrow = sapply(mlista,nrow)
  clone_thres = round(well_filter_thres * mean(alpha_nrow[alpha_nrow != 0]))
  #rm(mlist)

  qc<-get_good_wells_new(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  if(verbose) {
    msg1 = paste("Clone threshold for QC:", clone_thres)
    message(msg1)
    msg2 = paste("Alpha wells passing QC:", sum(qc$a))
    message(msg2)
    msg4 = paste("Beta wells passing QC:", sum(qc$b))
    message(msg4)
    msg3 = paste("Alpha wells failing QC:", sum(!qc$a))
    message(msg3)
    msg5 = paste("Beta wells failing QC:", sum(!qc$b))
    message(msg5)
  }

  #wells_a = sapply(strsplit(names(mlista), "_"), function(x) x[[well_pos]])
  #wells_b = sapply(strsplit(names(mlistb), "_"), function(x) x[[well_pos]])
  if(write_extra_files) {
    stats_a = tibble(
      a_names=df_meta$file_alpha,
      well = df_meta$well,
      a_sum_counts=sapply(mlista,function(x) { if(nrow(x) > 0) { return(x[,sum(readCount),]) } else { return(0) } }),
      a_rows=sapply(mlista,nrow),
      qc_pass_a=qc$a
    )
    stats_b = tibble(
      b_names=df_meta$file_beta,
      well = df_meta$well,
      b_sum_counts=sapply(mlistb,function(x) { if(nrow(x) > 0) { return(x[,sum(readCount),]) } else { return(0) } }),
      b_rows=sapply(mlistb,nrow),
      qc_pass_b=qc$b
    )
    plate_stats = full_join(stats_a, stats_b, by="well")
    # plate_stats = data.table::data.table(
    #   a_names=names(mlista),b_names=names(mlistb),
    #   a_sum_counts=sapply(mlista,function(x)x[,sum(readCount),]),
    #   b_sum_counts=sapply(mlistb,function(x)x[,sum(readCount),]),
    #   a_rows=sapply(mlista,nrow),
    #   b_rows=sapply(mlistb,nrow),
    #   qc_pass_a=qc$a,qc_pass_b=qc$b)
    fwrite(plate_stats,file.path(folder_out, paste0(prefix,"_plate_stats.tsv")),sep="\t")
  }

  # names(mlista) = df_meta$well
  # names(mlistb) = df_meta$well
  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  ## reorder wells
  # all_wells = get_well_subset(1:16,1:24)
  # sub_wells = all_wells[all_wells %in% qc$well_ids]
  # mlista = mlista[sub_wells]
  # mlistb = mlistb[sub_wells]

  if(exclude_nonfunctional) {
    message("Excluding chains with stop codons or frameshifts...")
    mlista = lapply(mlista, .get_functional)
    mlistb = lapply(mlistb, .get_functional)

    # tmpa2 = bind_rows(mlista)
    # tmpb2 = bind_rows(mlista)
    #
    # n_unique_alpha2 = length(unique(tmpa2$targetSequences))
    # n_unique_beta2 = length(unique(tmpb2$targetSequences))
    #
    # msg1 = paste(n_unique_alpha-n_unique_alpha2, "non-functional TCRalpha chains excluded")
    # msg2 = paste(n_unique_beta-n_unique_beta2, "non-functional TCRbeta chains excluded")
    # if(verbose) message(msg1)
    # if(verbose) message(msg2)

  }

  if(verbose) message("Tabulating TCRalpha pseudobulk counts")
  combd_a<-.combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  combd_a = combd_a[order(-readCount),]

  file_tra = paste0(prefix,"_pseudobulk_TRA.tsv")
  msg = paste("Writing TCRalpha pseudobulk file...", file_tra)
  if(verbose) message(msg)
  fwrite(combd_a,file.path(folder_out, file_tra),sep="\t")
  if(verbose) message("Tabulating TCRbeta pseudobulk counts")
  combd_b<-.combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  combd_b = combd_b[order(-readCount),]
  file_trb = paste0(prefix,"_pseudobulk_TRB.tsv")
  msg = paste("Writing TCRbeta pseudobulk file...", file_trb)
  if(verbose) message(msg)
  fwrite(combd_b,file.path(folder_out, file_trb),sep="\t")
  if(verbose) message("Pseudobulk done")

  if(verbose) message("Merging alpha clonesets...")
  bigma<-big_merge_freqs2(mlista, min_reads = min_reads)
  if(verbose) message(paste("Done! Unique alpha clones after filtering:", nrow(bigma)))
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  if(verbose) message(paste0("Unique alpha clones in more than ",min_wells," wells: ", nrow(bigmas)))
  bigmb<-big_merge_freqs2(mlistb, min_reads = min_reads)
  if(verbose) message(paste("Done! Unique beta clones after filtering:", nrow(bigmb)))
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  if(verbose) message(paste0("Unique beta clones and wells in more than ",min_wells," wells: ", nrow(bigmbs)))
  if(write_extra_files) {
    if(verbose) message("Writing TCRalpha/beta read count matrices (well x chain) and single-chain metadata...")
    check = all.equal(colnames(bigmas), colnames(bigmbs))
    if(!check) stop("Error creating alpha/beta read count matrices: wells are not in the same order")
    alpha_meta = combd_a[match(rownames(bigmas), combd_a$targetSequences),]
    beta_meta = combd_b[match(rownames(bigmbs), combd_b$targetSequences),]
    col_data = data.frame(idx = 1:length(colnames(bigmas)),well = colnames(bigmas))

    file_alpha_mat = paste(prefix,"alpha_mat.rds",sep="_")
    file_beta_mat = paste(prefix,"beta_mat.rds",sep="_")
    file_alpha_meta = paste(prefix,"alpha_meta.parquet",sep="_")
    file_beta_meta = paste(prefix,"beta_meta.parquet",sep="_")
    file_well_meta = paste(prefix,"well_meta.parquet",sep="_")

    if(verbose) message(paste("Writing TCRalpha read count matrix...", file_alpha_mat))
    saveRDS(t(bigmas), file = file.path(folder_out, file_alpha_mat)) ## columns are clones
    if(verbose) message(paste("Writing TCRbeta read count matrix...", file_beta_mat))
    saveRDS(t(bigmbs), file = file.path(folder_out, file_beta_mat)) ## columns are clones
    if(verbose) message(paste("Writing chain metadata for TCRalpha matrix...", file_alpha_meta))
    nanoparquet::write_parquet(alpha_meta,file.path(folder_out, file_alpha_meta))
    if(verbose) message(paste("Writing chain metadata for TCRbeta matrix...", file_beta_meta))
    nanoparquet::write_parquet(beta_meta,file.path(folder_out, file_beta_meta))
    if(verbose) message(paste("Writing well metadata for read count matrices...", file_well_meta))
    nanoparquet::write_parquet(col_data,file.path(folder_out, file_well_meta))
    #write(rownames(bigmas),file=file.path(folder_out, paste0(prefix,"_bigmas_names.tsv")))
    #write(rownames(bigmbs),file=file.path(folder_out, paste0(prefix,"_bigmbs_names.tsv")))
    # write_dat(as.matrix(bigmas),fname = file.path(folder_out, paste0(prefix,"_bigmas.tsv")))
    # write_dat(as.matrix(bigmbs),fname = file.path(folder_out, paste0(prefix,"_bigmbs.tsv")))
    # write.table(summary(bigmas), file.path(folder_out, paste0(prefix,"_bigmas_sparse_coo.tsv")), sep = "\t", row.names = FALSE)
    # write.table(summary(bigmbs), file.path(folder_out, paste0(prefix,"_bigmbs_sparse_coo.tsv")), sep = "\t", row.names = FALSE)
  }
  n_wells=ncol(bigmas)
  if(verbose) message("Pre-computing look-up table:")
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  if(write_extra_files) write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))

  if(compute==T) {
    message("Running pairing algorithms...")
    # reticulate stuff =======================
    py_path = system.file("python/pairing/", package = "TIRTLtools")
    pairing = reticulate::import_from_path("pairing_all_backends", path = py_path, convert = TRUE, delay_load = TRUE)

    bigmas_py = np_array(as.matrix(bigmas), dtype = "float32")
    bigmbs_py = np_array(as.matrix(bigmbs), dtype = "float32")
    mdh_py = r_to_py(mdh)

    pair_res = pairing$pairing(
      prefix = prefix, folder_out = folder_out, bigmas = bigmas_py, bigmbs = bigmbs_py,
      mdh = mdh_py, backend = backend, filter_before_top3 = filter_before_top3, chunk_size = chunk_size,
      write_files = write_extra_files
      )
  }
  message("Pairing is finished.")
  ### fix for when r_to_py doesn't convert data frames
  if(reticulate::is_py_object(pair_res[[1]])) pair_res = .fix_py_to_r_df_list(pair_res)

  print("Filtering results, adding amino acid and V segment information...")

  n_wells = ncol(bigmas)

  gpu_res = pair_res$mdh %>%
    mutate(alpha_nuc_seq = rownames(bigmas)[alpha_nuc],
           beta_nuc_seq = rownames(bigmbs)[beta_nuc]
    ) %>%
    mutate(alpha_nuc = alpha_nuc_seq, beta_nuc = beta_nuc_seq) %>%
    mutate(alpha_beta = paste0(alpha_nuc_seq,"_",beta_nuc_seq)) %>%
    mutate(method = "madhype") %>%
    as.data.table()
  gpu_res_corr = pair_res$corr %>%
    mutate(alpha_nuc_seq = rownames(bigmas)[alpha_nuc],
           beta_nuc_seq = rownames(bigmbs)[beta_nuc]) %>%
    mutate(alpha_nuc = alpha_nuc_seq, beta_nuc = beta_nuc_seq) %>%
    mutate(alpha_beta = paste0(alpha_nuc_seq,"_",beta_nuc_seq)) %>%
    mutate(ts = r*sqrt((n_wells - 2) / (1 - r^2))) %>%
    mutate(pval = 2*pt(-abs(ts), n_wells - 2)) %>%
    as.data.table()

  gpu_res_corr[,pval_adj:=pval/sort(pval)[3],alpha_nuc]
  gpu_res_corr[,method:="tshell",]

  # res_gpu$ts=res_gpu$r*sqrt((n_wells - 2) / (1 - res_gpu$r^2))
  # res_gpu$pval=2 * pt(-abs(res_gpu$ts), n_wells - 2)
  # res_gpu[,pval_adj:=pval/sort(pval)[3],alpha_nuc]
  # res_gpu[,method:="tshell",]

  # I also want to compute
  result<-rbind(gpu_res,gpu_res_corr,fill=T)
  #probably this we don't want. what is it???
  #result<-as.data.table(result[!duplicated(result[,1:2]),])

  result[,loss_a_frac:=(wb-wij)/(wij+(wb-wij)+(wa-wij)),]
  result[,loss_b_frac:=(wa-wij)/(wij+(wb-wij)+(wa-wij)),]
  result[,wi:=wa-wij,]
  result[,wj:=wb-wij,]

  #groom the output.

  unique_combinations <- unique(result[, .(wi, wj, wij)])
  print("Scoring unique pairs...")
  for (i in 1:nrow(unique_combinations)){
    #if(i%%1000==0)print(i)
    unique_combinations$score[i]<-log10(estimate_pair_prob(wi = unique_combinations[i,]$wi,wj = unique_combinations[i,]$wj,w_ij = unique_combinations[i,]$wij,n_wells,cpw = clone_thres,alpha = 2,prior=1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5))
  }

  #result<-result[order(-method),][!duplicated(alpha_beta),]
  #result<-merge(result, unique_combinations, by = c("wi", "wj", "wij"), all.x = TRUE)[method=="madhype"|(`method`=="tshell"&`wij`>get("wij_thres_tshell")&`pval_adj`<get("pval_thres_tshell")&(`loss_a_frac`+`loss_b_frac`)<0.5),]#there was no filter here before, check if it works!!!
  merged <- merge(result, unique_combinations, by = c("wi", "wj", "wij"), all.x = TRUE)
  result <- merged[method=="madhype"|(`method`=="tshell"&`wij`>wij_thres_tshell&`pval_adj`<pval_thres_tshell&(`loss_a_frac`+`loss_b_frac`)<0.5),]

  #result<-result[(((loss_a_frac+loss_b_frac)<0.5)&(wij>3))|(score>0.1),]

  tp_a<-add_VJ_aa(result$alpha_nuc,rbindlist(mlista))
  result$cdr3a=tp_a$cdr3aa
  result$va=tp_a$v
  result$ja=tp_a$j

  tp_b<-add_VJ_aa(result$beta_nuc,rbindlist(mlistb))
  result$cdr3b=tp_b$cdr3aa
  result$vb=tp_b$v
  result$jb=tp_b$j

  ### (optional) secondary algorithm to select best MAD-HYPE/T-SHELL clones
  result = .filter_pairs(result, filter_madhype = select_best_madhype, filter_tshell = select_best_tshell)

  file_paired = paste0(prefix,"_TIRTLoutput.tsv")
  if(verbose) message(paste("Writing TCRalpha/beta pairs...", file_paired))
  fwrite(result, file.path(folder_out, file_paired),sep="\t")
  if(verbose) {
    message("All pairing is finished!")
    message(paste("Number of clones paired by MAD-HYPE algorithm:", sum(result$method == "madhype")))
    message(paste("Number of clones paired by T-SHELL algorithm:", sum(result$method == "tshell")))
    message(paste("Total number of unique TCRalpha/beta pairs:", length(unique(result$alpha_beta))))
  }
  tictoc::toc()
  return(result)
}

.combineTCR<-function(dt){
  dt[, v := tstrsplit(allVHitsWithScore, "*", fixed = TRUE, fill = "")[[1]]]
  dt[, j := tstrsplit(allJHitsWithScore, "*", fixed = TRUE, fill = "")[[1]]]
  setkey(dt, targetSequences)
  nmax=length(unique(dt$file))
  #print(nmax)
  out <- dt[,
            .(
              readCount = sum(readCount),
              v = v[which.max(tabulate(match(v, v)))],
              j = j[which.max(tabulate(match(j, j)))],
              aaSeqCDR3 = aaSeqCDR3[which.max(tabulate(match(aaSeqCDR3, aaSeqCDR3)))],
              n_wells=.N,
              readCount_max=max(readCount),
              readCount_median=median(readCount),
              avg=sum(readFraction)/nmax,
              sem=sd(c(readFraction,rep(0,times=nmax-.N)))/sqrt(nmax)
            ),
            by = targetSequences #or by EACHI?
  ]
  out[, readFraction := readCount / sum(readCount)]
  out  #get most popular V and J and aaSeqCDR3 for each target
}

.get_functional = function(df) {
  check1 = !grepl("\\*", df$aaSeqCDR3)
  check2 = !grepl("_", df$aaSeqCDR3)
  df[check1 & check2,]
}

## sort by score (madhype) or p-value (tshell) and select best pairs
.filter_pairs = function(df, filter_madhype = TRUE, filter_tshell = TRUE) {
  df = df %>%
    mutate(is_functional = (!grepl("\\*", cdr3a)) & (!grepl("_", cdr3a)) & (!grepl("\\*", cdr3b)) & (!grepl("_", cdr3b)) )

  ### make final pairs data frame after sorting by score (or p-value for t-shell)
  .build_pairs_df = function(df) {
    pairs = df[c(),]
    for(i in 1:nrow(df)) {
      #print(i)
      alpha_tmp = df$alpha_nuc[i]
      beta_tmp = df$beta_nuc[i]
      n_beta_paired = sum(pairs$beta_nuc %in% beta_tmp)
      n_alpha_paired = sum(pairs$alpha_nuc %in% alpha_tmp)
      if(n_beta_paired <= 1 && n_alpha_paired == 0) {
        tmp = df[i,]
        pairs = bind_rows(pairs, tmp)
      } ## else discard row
    }
    pairs = pairs %>% select(-is_functional)
    return(pairs)
  }

  df_mdh = df %>% filter(method == "madhype") %>%
    arrange(desc(is_functional), desc(score))
  df_tshell = df %>% filter(method == "tshell") %>%
    arrange(desc(is_functional), pval)
  if(filter_madhype) {
    pairs_mdh = .build_pairs_df(df_mdh)
  } else {
    pairs_mdh = df_mdh
  }
  if(filter_tshell) {
    pairs_tshell = .build_pairs_df(df_tshell)
  } else {
    pairs_tshell = df_tshell
  }
  pairs_both = bind_rows(pairs_mdh, pairs_tshell)
  return(pairs_both)
}
