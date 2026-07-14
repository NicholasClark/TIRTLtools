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
#' @param tshell_settings a named list with pairing settings for T-SHELL: `pval_thres_tshell` and `wij_thres_tshell`. You may pass your own list or call \code{\link{get_tshell_settings}()} with setting "auto", "384_well", or "96_well". By default, "auto" is selected and "384_well" settings will be used if there are >= 150 wells passing QC, otherwise "96_well" settings will be used.
#' - `pval_thres_tshell` is the adjusted p-value threshold for T-SHELL significance (default 1e-10 for 384-well plate, 1e-3 for 96-well plate)
#' - `wij_thres_tshell` is the threshold for the number of wells containing both chains for T-SHELL significance (default >2 wells for 384-well plate, >3 wells for 96-well plate)
#' @param wellset a vector of wells to use for the pairing
#' @param backend the computing backend to use. The function looks for a GPU and automatically chooses an appropriate backend by default.
#' @param verbose whether to print out messages (default TRUE)
#' @param write_extra_files whether to write un-necessary intermediate files (default FALSE)
#' @param filter_before_top3 whether to filter by loss fraction before extracting top 3 correlation values for T-SHELL (default FALSE)
#' @param chunk_size batch size for calculations in pairing scripts
#' @param exclude_nonfunctional whether to exclude non-functional chains before pairing (default is FALSE)
#' @param select_best_madhype whether to use a secondary algorithm on the pairs from the MAD-HYPE algorithm to select
#' the best pairs for each clone (default is FALSE)
#' @param select_best_tshell whether to use a secondary algorithm on the pairs from the T-SHELL algorithm to select
#' the best pairs for each clone (default is FALSE)
#' @param gzip_output whether to compress (gzip) output files (default is FALSE)
#' @param pseudobulk_only only output pseudobulk with no pairing (default is FALSE)
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
    tshell_settings = get_tshell_settings(format = "auto"),
    wellset = get_well_subset(1:16,1:24),
    wellset1 = lifecycle::deprecated(),
    compute = lifecycle::deprecated(),
    backend = c("auto", "cpu", "cupy", "mlx"),
    pval_thres_tshell = lifecycle::deprecated(),
    wij_thres_tshell= lifecycle::deprecated(),
    verbose = TRUE,
    write_extra_files = FALSE,
    filter_before_top3 = FALSE,
    fork = lifecycle::deprecated(),
    shared = lifecycle::deprecated(),
    chunk_size = 500,
    exclude_nonfunctional = FALSE,
    select_best_madhype = FALSE,
    select_best_tshell = FALSE,
    gzip_output = FALSE,
    pseudobulk_only = FALSE
){
  tictoc::tic()

  if (lifecycle::is_present(wellset1)) {
    lifecycle::deprecate_warn(
      when = "0.1.22",
      what = "run_pairing(wellset1=)",
      with = "run_pairing(wellset=)"
    )
    wellset = wellset1
  }
  if(!is.null(tshell_settings)) {
    checkmate::assert_subset(c("pval_thres_tshell", "wij_thres_tshell"), names(tshell_settings))
    checkmate::assert_number(tshell_settings$pval_thres_tshell, lower = 0, upper = 1)
    checkmate::assert_int(tshell_settings$wij_thres_tshell, lower = 0, upper = 383)
  }

  #ensure_python_env()
  py_require( packages = get_py_deps() )

  backend = backend[1]

  .create_folder(folder_out)

  df_alpha = tibble(file_alpha = list.files(path = folder_path,full.names = F, pattern = ".*_TRA\\.tsv")) %>%
    mutate(file_alpha_no_dot = gsub(".","_",file_alpha,fixed=T),
           #well = strsplit(file_alpha_no_dot, "_") %>% sapply(., function(x) x[[well_pos]]),
           well = .get_well_from_filename(file_alpha_no_dot),
           well_row = substr(well, 0, 1),
           well_column = substr(well,2,3) %>% as.integer(),
           in_wellset = well %in% wellset) %>%
    filter(in_wellset) %>%
    arrange(well_row, well_column) %>%
    select(-well_row, -well_column, -in_wellset)
  df_beta = tibble(file_beta = list.files(path = folder_path,full.names = F, pattern = ".*_TRB\\.tsv")) %>%
    mutate(file_beta_no_dot = gsub(".","_",file_beta,fixed=T),
           #well = strsplit(file_beta_no_dot, "_") %>% sapply(., function(x) x[[well_pos]]),
           well = .get_well_from_filename(file_beta_no_dot),
           well_row = substr(well, 0, 1),
           well_column = substr(well,2,3) %>% as.integer(),
           in_wellset = well %in% wellset)  %>%
    filter(in_wellset) %>%
    arrange(well_row, well_column) %>%
    select(-well_row, -well_column, -in_wellset)
  missing_alpha_wells = wellset[!wellset %in% df_alpha$well]
  missing_beta_wells = wellset[!wellset %in% df_beta$well]
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
  wells_missing = wellset[!wellset %in% df_meta$well]
  df_meta = df_meta %>%
    bind_rows(tibble(well = wells_missing)) %>%
    mutate(well_row = substr(well, 0, 1), well_column = substr(well,2,3) %>% as.integer()) %>%
    arrange(well_row, well_column)
  # missing_wells_final = wellset[!wellset %in% df_meta$well]
  # if(length(missing_wells_final) > 0) {
  #   msg = paste(length(missing_wells_final), "wells excluded for missing .tsv files:", paste0(missing_wells_final, collapse = ", "))
  #   if(verbose) message(msg)
  # }

  # files = list.files(path = folder_path,full.names = F)
  # files_no_dot = gsub(".","_",files,fixed=T)
  # file_wells = strsplit(files_no_dot, "_") %>% sapply(., function(x) x[[well_pos]])
  # files_bool = file_wells %in% wellset
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
  beta_nrow = sapply(mlistb,nrow)
  clone_thres_alpha = round(well_filter_thres * mean(alpha_nrow[alpha_nrow != 0]))
  clone_thres_beta = round(well_filter_thres * mean(beta_nrow[beta_nrow != 0]))

  qc = get_good_wells_new(mlista,mlistb,
                         thres_alpha = clone_thres_alpha,
                         thres_beta = clone_thres_beta,
                         pos=well_pos,
                         wellset=wellset,
                         verbose = verbose)

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

  n_wells_pass = length(mlista)
  has_tshell_settings = sum(c("pval_thres_tshell", "wij_thres_tshell") %in% names(tshell_settings)) == 2
  if(!has_tshell_settings) {
    msg = "Using 'auto' T-SHELL settings, selecting 'pval_thres_tshell' and 'wij_thres_tshell' based on the number of wells passing QC"
    if(verbose) message(msg)
    n_well_thres = 150
    if(n_wells_pass >= n_well_thres) {
      tshell_settings = get_tshell_settings("384_well")
      msg = paste(n_wells_pass, "wells passing: Using '384_well' T-SHELL settings")
      if(verbose) message(msg)
    } else {
      tshell_settings = get_tshell_settings("96_well")
      msg = paste(n_wells_pass, "wells passing: Using '96_well' T-SHELL settings")
      if(verbose) message(msg)
    }
  }

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
  if(gzip_output) file_tra = paste0(file_tra, ".gz")
  msg = paste("Writing TCRalpha pseudobulk file...", file_tra)
  if(verbose) message(msg)
  fwrite(combd_a,file.path(folder_out, file_tra),sep="\t")
  if(verbose) message("Tabulating TCRbeta pseudobulk counts")
  combd_b<-.combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  combd_b = combd_b[order(-readCount),]
  file_trb = paste0(prefix,"_pseudobulk_TRB.tsv")
  if(gzip_output) file_trb = paste0(file_trb, ".gz")
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

  if(!pseudobulk_only) {
    if(verbose) message("Pre-computing look-up table:")
    mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres_alpha,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
    if(write_extra_files) write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))

    #if(compute==T) {
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
        ) ## returns a list with two data frames: "mdh" (mad-hype) and "corr" (T-shell)
    #}
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
    ## T-shell results -- top 3 betas for each alpha
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
    result <- merged[method=="madhype"|(`method`=="tshell"&`wij`>tshell_settings$wij_thres_tshell&`pval_adj`<tshell_settings$pval_thres_tshell&(`loss_a_frac`+`loss_b_frac`)<0.5),]

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
    if(gzip_output) file_paired = paste0(file_paired, ".gz")
    fwrite(result, file.path(folder_out, file_paired),sep="\t")
    if(verbose) {
      message("All pairing is finished!")
      message(paste("Number of clones paired by MAD-HYPE algorithm:", sum(result$method == "madhype")))
      message(paste("Number of clones paired by T-SHELL algorithm:", sum(result$method == "tshell")))
      message(paste("Total number of unique TCRalpha/beta pairs:", length(unique(result$alpha_beta))))
    }
  } else {
    result = NULL
  }
  tictoc::toc()
  if(is.null(result)) return(invisible(NULL))
  return(result)
}
