run_pairing = function(
    folder_path,
    folder_out,
    prefix="tmp",
    well_filter_thres=0.75,
    min_reads=0,
    min_wells=2,
    well_pos=3,
    wellset1=get_well_subset(1:16,1:24),
    compute=T,
    backend = c("auto", "cpu", "cupy", "mlx"),
    pval_thres_tshell=1e-10,
    wij_thres_tshell=2,
    verbose = TRUE,
    write_extra_files = FALSE,
    filter_before_top3 = FALSE,
    fork = NULL,
    shared = NULL
){
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

  mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  if(verbose) {
    print("clonesets loaded")
    print(names(mlist))
    print(Sys.time())
  }
  wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  if(verbose) {
    print("clone_threshold")
    print(clone_thres)
    print("alpha_wells_working")
    print(table(qc$a))
    print("beta_wells_working")
    print(table(qc$b))
  }

  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  if(verbose) {
    print("pseudobulk_alpha")
    print(Sys.time())
  }
  combd_a<-.combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  if(verbose) print("pseudobulk_beta")
  fwrite(combd_a[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRA.tsv")),sep="\t")
  if(verbose) print(Sys.time())
  combd_b<-.combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRB.tsv")),sep="\t")
  if(verbose) print("pseudobulk_done")

  if(verbose) {
    print(Sys.time())
    print("start")
  }
  bigma<-big_merge_freqs2(mlista,min_reads = min_reads)
  if(verbose) print(dim(bigma))
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  if(verbose) {
    print(Sys.time())
    print("big merge done")
    print("bigmas")
    print(dim(bigmas))
  }
  bigmb<-big_merge_freqs2(mlistb,min_reads = min_reads)
  if(verbose) print(dim(bigmb))
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  if(verbose) {
    print(Sys.time())
    print("big merge done")
    print("bigmbs")
    print(dim(bigmbs))
  }
  if(write_extra_files) {
    write(rownames(bigmas),file=file.path(folder_out, paste0(prefix,"_bigmas_names.tsv")))
    write(rownames(bigmbs),file=file.path(folder_out, paste0(prefix,"_bigmbs_names.tsv")))
    write_dat(as.matrix(bigmas),fname = file.path(folder_out, paste0(prefix,"_bigmas.tsv")))
    write_dat(as.matrix(bigmbs),fname = file.path(folder_out, paste0(prefix,"_bigmbs.tsv")))
  }
  if(verbose) print(Sys.time())
  n_wells=ncol(bigmas)
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  if(write_extra_files) write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))
  if(verbose) print(Sys.time())

  if(compute==T) {
    #### basilisk stuff
    is_loaded = pkgload::is_dev_package("TIRTLtools")
    if(is_loaded) {
      if(is.null(fork)) fork = FALSE
      if(is.null(shared)) shared = TRUE
      print("package loaded by pkgload/devtools")
    } else {
      if(is.null(fork)) fork = FALSE
      if(is.null(shared)) shared = FALSE
      print("package loaded from installed version")
    }
    #env = .choose_basilisk_env(backend)
    env = TIRTLtools_py_env
    proc <- basilisk::basiliskStart(env, fork = fork, shared = shared)
    on.exit(basilisk::basiliskStop(proc))
    py_path = system.file("python/pairing/", package = "TIRTLtools")

    pair_res = basilisk::basiliskRun(proc, fun=function(prefix, folder_out, bigmas, bigmbs, mdh, backend, filter_before_top3) {
      pairing = reticulate::import_from_path("pairing_all_backends", path = py_path, convert = TRUE, delay_load = TRUE)

      bigmas_py = np_array(as.matrix(bigmas), dtype = "float32")
      bigmbs_py = np_array(as.matrix(bigmbs), dtype = "float32")
      mdh_py = r_to_py(mdh)

      pairing_res = pairing$pairing(prefix = prefix, folder_out = folder_out,
                                 bigmas = bigmas_py, bigmbs = bigmbs_py, mdh = mdh_py,
                                 backend = backend, filter_before_top3 = filter_before_top3)
      return(pairing_res)
    }, prefix = prefix, folder_out = folder_out, bigmas = bigmas, bigmbs = bigmbs, mdh = mdh, backend = backend, filter_before_top3 = filter_before_top3)

  }

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
  for (i in 1:nrow(unique_combinations)){
    if(i%%1000==0)print(i)
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
