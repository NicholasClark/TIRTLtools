

run_single_point_analysis_sub_gpu_pseudobulk = function(
    folder_path,
    folder_out,
    prefix="tmp",
    well_filter_thres=0.75,
    min_reads=0,
    min_wells=2,
    well_pos=3,
    wellset1=get_well_subset(1:16,1:24),
    compute=T,
    pval_thres_tshell=1e-10,
    wij_thres_tshell=2,
    verbose = TRUE){

  if(verbose) {
    message("start")
    message(Sys.time())
  }

  ## Create folder for results
  if (!dir.exists(folder_out)) {
    dir.create(folder_out, recursive = TRUE)
    msg = paste("Created folder:", folder_out)
    if(verbose) message(msg)
  } else {
    msg = paste("Folder already exists:", folder_out)
    if(verbose) message(msg)
  }

  mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  if(verbose) {
    message("clonesets loaded")
    message(names(mlist))
    message(Sys.time())
  }
  wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  if(verbose) {
    message("clone_threshold")
    message(clone_thres)
    message("alpha_wells_working")
    message(table(qc$a))
    message("alpha_wells_working")
    message(table(qc$a))
  }

  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  message("pseudobulk_alpha")
  message(Sys.time())
  combd_a<-combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  message("pseudobulk_beta")
  fwrite(combd_a[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRA.tsv"),sep="\t"))
  message(Sys.time())
  combd_b<-combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b[order(-readCount),],file.path(folder_out, paste0(prefix,"_pseudobulk_TRB.tsv"),sep="\t"))
  message("pseudobulk_done")

  message(Sys.time())
  message("start")
  bigma<-big_merge_freqs2(mlista,min_reads = min_reads)
  message(dim(bigma))
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  message(Sys.time())
  message("big merge done")
  message("bigmas")
  message(dim(bigmas))
  bigmb<-big_merge_freqs2(mlistb,min_reads = min_reads)
  message(dim(bigmb))
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  message(Sys.time())
  message("big merge done")
  message("bigmbs")
  message(dim(bigmbs))
  write(rownames(bigmas),file=file.path(folder_out, paste0(prefix,"_bigmas_names.tsv")))
  write(rownames(bigmbs),file=file.path(folder_out, paste0(prefix,"_bigmbs_names.tsv")))
  write_dat(as.matrix(bigmas),fname = file.path(folder_out, paste0(prefix,"_bigmas.tsv")))
  write_dat(as.matrix(bigmbs),fname = file.path(folder_out, paste0(prefix,"_bigmbs.tsv")))
  message(Sys.time())
  n_wells=ncol(bigmas)
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))
  message(Sys.time())



  if(compute==T)system(paste0("python3 cupy_madhype_script.py ",prefix,collapse=""))
  # here goes SYS call to python script.
  #python3 mlx_madhype_script ~/R_projects/mlx_dev/plate6
  # and here we go read it:
  gpu_res<-read_gpu(prefix)
  gpu_res_corr<-read_gpu_corr(prefix)
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
    if(i%%1000==0)message(i)
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

  fwrite(result, file.path(folder_out, paste0(prefix,"_TIRTLoutput.tsv"),sep="\t"))
  return(result)
}

run_single_point_analysis_sub_gpu_pseudobulk_legacy<-function(folder_path,prefix="tmp",well_filter_thres=0.75,min_reads=0,min_wells=2,well_pos=3,wellset1=get_well_subset(1:16,1:24),compute=T,pval_thres_tshell=1e-10,wij_thres_tshell=2){ #this is with gpu backend
  print("start")
  print(Sys.time())
  mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  print("clonesets loaded")
  print(names(mlist))
  print(Sys.time())
  wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  print("clone_threshold")
  print(clone_thres)
  print("alpha_wells_working")
  print(table(qc$a))
  print("alpha_wells_working")
  print(table(qc$a))
  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  print("pseudobulk_alpha")
  print(Sys.time())
  combd_a<-combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  print("pseudobulk_beta")
  fwrite(combd_a[order(-readCount),],paste0(prefix,"_pseudobulk_TRA.tsv"),sep="\t")
  print(Sys.time())
  combd_b<-combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b[order(-readCount),],paste0(prefix,"_pseudobulk_TRB.tsv"),sep="\t")
  print("pseudobulk_done")

    print(Sys.time())
  print("start")
  bigma<-big_merge_freqs2(mlista,min_reads = min_reads)
  print(dim(bigma))
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  print(Sys.time())
  print("big merge done")
  print("bigmas")
  print(dim(bigmas))
  bigmb<-big_merge_freqs2(mlistb,min_reads = min_reads)
  print(dim(bigmb))
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  print(Sys.time())
  print("big merge done")
  print("bigmbs")
  print(dim(bigmbs))
  write(rownames(bigmas),file=paste0(prefix,"_bigmas_names.tsv"))
  write(rownames(bigmbs),file=paste0(prefix,"_bigmbs_names.tsv"))
  write_dat(as.matrix(bigmas),fname = paste0(prefix,"_bigmas.tsv"))
  write_dat(as.matrix(bigmbs),fname = paste0(prefix,"_bigmbs.tsv"))
  print(Sys.time())
  n_wells=ncol(bigmas)
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  write_dat(mdh,fname = paste0(prefix,"_mdh.tsv"))
  print(Sys.time())
  if(compute==T)system(paste0("python3 cupy_madhype_script.py ",prefix,collapse=""))
  # here goes SYS call to python script.
  #python3 mlx_madhype_script ~/R_projects/mlx_dev/plate6
  # and here we go read it:
  gpu_res<-read_gpu(prefix)
  gpu_res_corr<-read_gpu_corr(prefix)
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

  fwrite(result,paste0(prefix,"_TIRTLoutput.tsv"),sep="\t")
  return(result)
}

run_single_point_analysis_sub_gpu_pseudobulk_only2<-function(folder_path,prefix="tmp",well_filter_thres=0.75,min_reads=0,min_wells=2,well_pos=3,wellset1=get_well_subset(1:16,1:24),compute=T,pval_thres_tshell=1e-10,wij_thres_tshell=2){ #this is with gpu backend
  print("start")
  print(Sys.time())
  fullnames<-list.files(path = folder_path,full.names = T)
  mlist<-lapply(fullnames,fread)
  names(mlist)<-gsub(".","_",basename(fullnames),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  print("clonesets loaded")
  print(names(mlist))
  print(Sys.time())
  wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  print("clone_threshold")
  print(clone_thres)
  print("alpha_wells_working")
  print(table(qc$a))
  print("alpha_wells_working")
  print(table(qc$a))
  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  print(names(mlista[qc$a]))

  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc
  print("pseudobulk_alpha")
  print(Sys.time())
  combd_a<-combineTCR(rbindlist(mlista,idcol="file"))
  combd_a$max_wells<-sum(qc$a)
  print("pseudobulk_beta")
  fwrite(combd_a[order(-readCount),],paste0(prefix,"_pseudobulk_TRA.tsv"),sep="\t")
  print(Sys.time())
  combd_b<-combineTCR(rbindlist(mlistb,idcol="file"))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b[order(-readCount),],paste0(prefix,"_pseudobulk_TRB.tsv"),sep="\t")
  print("pseudobulk_done")
}



run_single_point_analysis_sub_pseudobulk_only<-function(folder_path,prefix="tmp",well_filter_thres=0.5,min_reads=0,min_wells=2,well_pos=3,wellset1=get_well_subset(1:16,1:24)){ #this is with gpu backend
  print("start")
  print(Sys.time())
  mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  print("clonesets loaded")
  print(Sys.time())
  wellsub<-sapply(strsplit(names(mlist),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlist,nrow)[wellsub]))
  rm(mlist)

  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  print("clone_threshold")
  print(clone_thres)
  print("alpha_wells_working")
  print(table(qc$a))
  print("alpha_wells_working")
  print(table(qc$a))
  #result<-do_analysis_madhyper_r_optim_both(mlista[qc$a],mlistb[qc$b],n_cells = clone_thres)
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc

  print("pseudobulk_alpha")
  print(Sys.time())
  combd_a<-combineTCR(rbindlist(mlista[qc$a])) #this is a mistake!
  combd_a$max_wells<-sum(qc$a)
  print("pseudobulk_beta")
  fwrite(combd_a,paste0(prefix,"_pseudobulk_TRA.tsv"),sep="\t")
  print(Sys.time())
  combd_b<-combineTCR(rbindlist(mlistb[qc$b]))
  combd_b$max_wells<-sum(qc$b)
  fwrite(combd_b,paste0(prefix,"_pseudobulk_TRB.tsv"),sep="\t")
  print("pseudobulk_done")
}

combineTCR<-function(dt){
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
