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
#' @param exclude_nonfunctional whether to exclude non-functional chains before pairing (default FALSE)
#'
#' @return
#' A data frame with the TCR-alpha/TCR-beta pairs.
#'
#' The function also writes three files to "folder_out":
#'  - A data frame ("<prefix>_pseudobulk_TRA.tsv") of pseudobulk counts for TCRalpha chains
#'  - A data frame ("<prefix>_pseudobulk_TRB.tsv") of pseudobulk counts for TCRbeta chains
#'  - A data frame ("<prefix>_TIRTLoutput.tsv") of TCR-alpha/TCR-beta pairs.
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
    exclude_nonfunctional = FALSE
){
  #ensure_python_env()
  py_require( packages = get_py_deps() )

  backend = backend[1]
  if(verbose) {
    print("start")
    print(Sys.time())
  }

  ## Create folder for results
  if (!dir.exists(folder_out)) {
    dir.create(folder_out, recursive = TRUE)
    msg = paste("Created folder:", folder_out)
    if(verbose) print(msg)
  } else {
    msg = paste("Folder already exists:", folder_out)
    if(verbose) print(msg)
  }

  files = list.files(path = folder_path,full.names = F)
  files_no_dot = gsub(".","_",files,fixed=T)
  file_wells = strsplit(files_no_dot, "_") %>% sapply(., function(x) x[[well_pos]])
  files_bool = file_wells %in% wellset1
  msg = paste("Reading", sum(files_bool), "wells from", folder_path)
  if(verbose) print(msg)
  files_load = file.path(folder_path, files[files_bool])
  mlist<-lapply(1:length(files_load),function(i){
    print(i)
    fread(files_load[i])
  })
  names(mlist) = files_no_dot[files_bool]

  #mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  #names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  if(verbose) {
    n_filesA = length(mlista)
    n_filesB = length(mlistb)
    msgA = paste(n_filesA, "TCRalpha well files loaded")
    msgB = paste(n_filesB, "TCRbeta well files loaded")
    print(msgA)
    print(msgB)
    #print("clonesets loaded")
    #print(names(mlist))
    print(Sys.time())
  }
  #wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  #clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  #wellsub<-sapply(strsplit(names(mlista),split="_",fixed=T),"[[",well_pos)%in%wellset1
  #clone_thres = round(well_filter_thres * mean(sapply(mlista,nrow)[wellsub]))
  clone_thres = round(well_filter_thres * mean(sapply(mlista,nrow)))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  if(verbose) {
    print("Clone threshold for QC:")
    print(clone_thres)
    print("Alpha wells passing QC:")
    print(table(qc$a))
    print("Beta wells passing QC:")
    print(table(qc$b))
  }

  if(write_extra_files) {
    plate_stats = data.table::data.table(
      a_names=names(mlista),b_names=names(mlistb),
      a_sum_counts=sapply(mlista,function(x)x[,sum(readCount),]),
      b_sum_counts=sapply(mlistb,function(x)x[,sum(readCount),]),
      a_rows=sapply(mlista,nrow),
      b_rows=sapply(mlistb,nrow),
      qc_pass_a=qc$a,qc_pass_b=qc$b)
    fwrite(plate_stats,file.path(folder_out, paste0(prefix,"_plate_stats.tsv")),sep="\t")
  }

  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  if(exclude_nonfunctional) {
    mlista = lapply(mlista, .get_functional)
    mlistb = lapply(mlistb, .get_functional)
  }

  if(verbose) {
    print("Tabulating TCRalpha pseudobulk counts")
    print(Sys.time())
  }
  combd_a<-.combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  fwrite(combd_a[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRA.tsv")),sep="\t")
  if(verbose) {
    print("Tabulating TCRbeta pseudobulk counts")
    print(Sys.time())
  }
  combd_b<-.combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRB.tsv")),sep="\t")
  if(verbose) print("Pseudobulk done.")

  if(verbose) {
    print(Sys.time())
    print("Merging alpha clonesets...")
  }
  bigma<-big_merge_freqs2(mlista,min_reads = min_reads)
  if(verbose) {
    print("Done! Unique alpha clones and wells after filtering:")
    print(dim(bigma))
  }
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  if(verbose) {
    print(Sys.time())
    print(paste0("Unique alpha clones and wells in more than: ",min_wells," wells"))
    print(dim(bigmas))
  }
  bigmb<-big_merge_freqs2(mlistb,min_reads = min_reads)
  if(verbose) {
    print("Done! beta clones and wells after filtering:")
    print(dim(bigmb))
  }
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  if(verbose) {
    print(Sys.time())
    print(paste0("Unique beta clones and wells in more than: ",min_wells," wells"))
    print(dim(bigmbs))
  }
  if(write_extra_files) {
    write(rownames(bigmas),file=file.path(folder_out, paste0(prefix,"_bigmas_names.tsv")))
    write(rownames(bigmbs),file=file.path(folder_out, paste0(prefix,"_bigmbs_names.tsv")))
    write_dat(as.matrix(bigmas),fname = file.path(folder_out, paste0(prefix,"_bigmas.tsv")))
    write_dat(as.matrix(bigmbs),fname = file.path(folder_out, paste0(prefix,"_bigmbs.tsv")))
    write.table(summary(bigmas), file.path(folder_out, paste0(prefix,"_bigmas_sparse_coo.tsv")), sep = "\t", row.names = FALSE)
    write.table(summary(bigmbs), file.path(folder_out, paste0(prefix,"_bigmbs_sparse_coo.tsv")), sep = "\t", row.names = FALSE)
  }
  if(verbose) print(Sys.time())
  n_wells=ncol(bigmas)
  if(verbose) print("Pre-computing look-up table:")
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  if(write_extra_files) write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))
  if(verbose) print(Sys.time())

  if(compute==T) {
    # basilisk stuff =====================
    # is_loaded = pkgload::is_dev_package("TIRTLtools")
    # if(is_loaded) {
    #   if(is.null(fork)) fork = FALSE
    #   if(is.null(shared)) shared = TRUE
    #   #print("package loaded by pkgload/devtools")
    # } else {
    #   if(is.null(fork)) fork = FALSE
    #   if(is.null(shared)) shared = FALSE
    #   #print("package loaded from installed version")
    # }
    # #env = .choose_basilisk_env(backend)
    # env = TIRTLtools_py_env
    # proc <- basilisk::basiliskStart(env, fork = fork, shared = shared)
    # on.exit(basilisk::basiliskStop(proc))
    # py_path = system.file("python/pairing/", package = "TIRTLtools")
    #
    # pair_res = basilisk::basiliskRun(proc, fun=function(prefix, folder_out, bigmas, bigmbs, mdh, backend, filter_before_top3) {
    #   pairing = reticulate::import_from_path("pairing_all_backends", path = py_path, convert = TRUE, delay_load = TRUE)
    #
    #   bigmas_py = np_array(as.matrix(bigmas), dtype = "float32")
    #   bigmbs_py = np_array(as.matrix(bigmbs), dtype = "float32")
    #   mdh_py = r_to_py(mdh)
    #
    #   pairing_res = pairing$pairing(prefix = prefix, folder_out = folder_out,
    #                              bigmas = bigmas_py, bigmbs = bigmbs_py, mdh = mdh_py,
    #                              backend = backend, filter_before_top3 = filter_before_top3)
    #   return(pairing_res)
    # }, prefix = prefix, folder_out = folder_out, bigmas = bigmas, bigmbs = bigmbs, mdh = mdh, backend = backend, filter_before_top3 = filter_before_top3)

    # reticulate stuff =======================
    #pair_res = basilisk::basiliskRun(proc, fun=function(prefix, folder_out, bigmas, bigmbs, mdh, backend, filter_before_top3) {
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
      #return(pairing_res)
    #}, prefix = prefix, folder_out = folder_out, bigmas = bigmas, bigmbs = bigmbs, mdh = mdh, backend = backend, filter_before_top3 = filter_before_top3)

  }
  ### fix for when r_to_py doesn't convert data frames
  if(reticulate::is_py_object(pair_res[[1]])) pair_res = .fix_py_to_r_df_list(pair_res)

  print("Filtering results, adding amino acid and V segment information")

  # cupy_madhype_script = reticulate::import_from_path("cupy_madhype_script", path = system.file("python/pairing/", package = "TIRTLtools"), convert = TRUE, delay_load = TRUE)
  # cupy_madhype_script$madhyper_process(prefix)
  # cupy_madhype_script$correlation_process(prefix)
  # if(compute==T)system(paste0("python3 cupy_madhype_script.py ",prefix,collapse=""))
  # here goes SYS call to python script.
  #python3 mlx_madhype_script ~/R_projects/mlx_dev/plate6
  # and here we go read it:

  #gpu_res<-read_gpu(file.path(folder_out, prefix))
  #gpu_res_corr<-read_gpu_corr(file.path(folder_out, prefix))

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

  if(verbose) {
    print(Sys.time())
    print("All is done! Number of paired clones:")
    print(table(result$method))
  }
  fwrite(result, file.path(folder_out, paste0(prefix,"_TIRTLoutput.tsv")),sep="\t")
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
