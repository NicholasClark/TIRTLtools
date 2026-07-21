

process_well_data = function(
    folder_path,
    folder_out,
    prefix,
    well_filter_thres=0.5,
    min_reads=0,
    min_wells=2,
    well_pos=3,
    wellset1=get_well_subset(1:16,1:24),
    verbose = TRUE,
    exclude_nonfunctional = FALSE
    ) {
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

  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  if(verbose) {
    n_filesA = length(mlista)
    n_filesB = length(mlistb)
    msgA = paste(n_filesA, "TCRalpha well files loaded")
    msgB = paste(n_filesB, "TCRbeta well files loaded")
    print(msgA)
    print(msgB)
    print(Sys.time())
  }
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
  wells_a = sapply(strsplit(names(mlista), "_"), function(x) x[[well_pos]])
  wells_b = sapply(strsplit(names(mlistb), "_"), function(x) x[[well_pos]])
  stats_a = tibble(
    a_names=names(mlista),
    well = wells_a,
    a_sum_counts=sapply(mlista,function(x)x[,sum(readCount),]),
    a_rows=sapply(mlista,nrow),
    qc_pass_a=qc$a
    )
  stats_b = tibble(
    b_names=names(mlistb),
    well = wells_b,
    b_sum_counts=sapply(mlistb,function(x)x[,sum(readCount),]),
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

  write(rownames(bigmas),file=file.path(folder_out, paste0(prefix,"_bigmas_names.tsv")))
  write(rownames(bigmbs),file=file.path(folder_out, paste0(prefix,"_bigmbs_names.tsv")))
  write_dat(as.matrix(bigmas),fname = file.path(folder_out, paste0(prefix,"_bigmas.tsv")))
  write_dat(as.matrix(bigmbs),fname = file.path(folder_out, paste0(prefix,"_bigmbs.tsv")))
  write.table(summary(bigmas), file.path(folder_out, paste0(prefix,"_bigmas_sparse_coo.tsv")), sep = "\t", row.names = FALSE)
  write.table(summary(bigmbs), file.path(folder_out, paste0(prefix,"_bigmbs_sparse_coo.tsv")), sep = "\t", row.names = FALSE)

  n_wells=ncol(bigmas)
  if(verbose) print("Pre-computing look-up table:")
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  write_dat(mdh,fname = file.path(folder_out, paste0(prefix,"_mdh.tsv")))
  if(verbose) print(Sys.time())
}
