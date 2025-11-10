big_merge_freqs2<-function(dtlist,min_reads=0){ #optimised
  all_data <- rbindlist(dtlist, idcol = "source")[readCount>min_reads,]
  all_data$source <- factor(all_data$source, levels=unique(all_data$source))
  cols <- levels(all_data$source)
  all_data$targetSequences <- factor(all_data$targetSequences,levels=unique(all_data$targetSequences))
  rows <- levels(all_data$targetSequences)
  stmpsp<-sparseMatrix(
    i = as.integer(all_data$targetSequences),    # row indices
    j = as.integer(all_data$source),  # column indices
    x = all_data$readFraction,            # values
    dimnames = list(rows, cols)    # dimension names
  )
  stmpsp
}

geta<-function(dtlist){
  dtlist[grepl("TRA",names(dtlist))]
}

getb<-function(dtlist){
  dtlist[grepl("TRB",names(dtlist))]  
}

geta_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRA",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

getb_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRB",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

get_well_subset<-function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F))
}

get_good_wells_sub<-function(alpha_list,beta_list,thres,pos=4,wellset=get_well_subset(1:16,1:24)){
  wellinds_a<-(sapply(strsplit(names(alpha_list),"_"),"[[",pos))  
  wellinds_b<-(sapply(strsplit(names(beta_list),"_"),"[[",pos))  
  good_wells<-intersect(wellinds_a[sapply(alpha_list,nrow)>thres],wellinds_b[sapply(beta_list,nrow)>thres])
  list(a=wellinds_a%in%intersect(good_wells,wellset),b=wellinds_b%in%intersect(good_wells,wellset),well_ids=intersect(good_wells,wellset))
}

madhyper_surface<-function(n_wells,cells=1000,alpha=2,prior=1){
  new_cube=matrix(0,nrow=n_wells+1,ncol=n_wells+1)
  i=1
  pb <- txtProgressBar(min = 0, max = n_wells, style = 3)
  for (wij in (n_wells):0){
    setTxtProgressBar(pb, n_wells-wij)
    
    ans<-rep(0,n_wells+1)
    if (log10(estimate_pair_prob(wi=i-1,wj=0,w_ij=wij,w_tot = n_wells,cpw=cells,alpha=alpha,prior=prior))>0.1){
      while(log10(estimate_pair_prob(wi=i-1,wj=0,w_ij=wij,w_tot = n_wells,cpw=cells,alpha=alpha,prior=prior))>0.1){i=i+1}
    } 
    else
    {while((i>1)&(log10(estimate_pair_prob(wi=i-1,wj=0,w_ij=wij,w_tot = n_wells,cpw=cells,alpha=alpha,prior=prior))<0.1)){i=i-1}}  
    z=1
    for (j in i:1)
    {
      while(log10(estimate_pair_prob(wi=j-1,wj=z-1,w_ij=wij,w_tot = n_wells,cpw=cells,alpha=alpha,prior=prior))>0.1){z=z+1}
      ans[j]=z-1
    }
    new_cube[wij+1,]<-ans
  }
  close(pb)
  new_cube
}

concordance<-function(tirtlseq1,tirtlseq2)
{
  table(tirtlseq1[beta_nuc%in%tirtlseq2$beta_nuc,alpha_beta,]%in%tirtlseq2$alpha_beta)
}

write_for_gpu<-function(mlista,mlistb,n_cells=3000,alpha=2,min_reads=0,min_wells=2,prefix="")#writes out bigmas,bigmbs and mdh files. 
{
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
  write(rownames(bigmas),file=paste0(prefix,"bigmas_names.tsv"))
  write(rownames(bigmbs),file=paste0(prefix,"bigmbs_names.tsv"))
  write_dat(as.matrix(bigmas),fname = paste0(prefix,"bigmas.tsv"))
  write_dat(as.matrix(bigmbs),fname = paste0(prefix,"bigmbs.tsv"))
  print(Sys.time())
  n_wells=ncol(bigmas)
  mdh<-madhyper_surface(n_wells = n_wells,cells = n_cells,alpha=alpha,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  write_dat(mdh,fname = paste0(prefix,"mdh.tsv"))
}

#to test on a server
read_gpu<-function(prefix){
  res_gpu<-fread(paste0(prefix,"_madhyperesults.csv"))
  res_gpu_b<-fread(paste0(prefix,"_bigmbs_names.tsv"),header=F)
  res_gpu_a<-fread(paste0(prefix,"_bigmas_names.tsv"),header=F)
  res_gpu[,alpha_nuc_seq:=res_gpu_a$V1[alpha_nuc],]
  res_gpu[,beta_nuc_seq:=res_gpu_b$V1[beta_nuc],]
  res_gpu[,alpha_nuc:=alpha_nuc_seq,]
  res_gpu[,beta_nuc:=beta_nuc_seq,]
  res_gpu[,alpha_beta:=paste0(alpha_nuc_seq,"_",beta_nuc_seq),]
  res_gpu[,method:="madhype",]
  return(res_gpu)
}

read_gpu_corr<-function(prefix){
  res_gpu<-fread(paste0(prefix,"_corresults.csv"))
  res_gpu_b<-fread(paste0(prefix,"_bigmbs_names.tsv"),header=F)
  res_gpu_a<-fread(paste0(prefix,"_bigmas_names.tsv"),header=F)
  n_wells<-ncol(fread(paste0(prefix,"_bigmas.tsv"),header=F,nrows = 1))
  res_gpu[,alpha_nuc_seq:=res_gpu_a$V1[alpha_nuc],]
  res_gpu[,beta_nuc_seq:=res_gpu_b$V1[beta_nuc],]
  res_gpu[,alpha_nuc:=alpha_nuc_seq,]
  res_gpu[,beta_nuc:=beta_nuc_seq,]
  
  res_gpu[,alpha_beta:=paste0(alpha_nuc_seq,"_",beta_nuc_seq),]
  res_gpu$ts=res_gpu$r*sqrt((n_wells - 2) / (1 - res_gpu$r^2))
  res_gpu$pval=2 * pt(-abs(res_gpu$ts), n_wells - 2)
  res_gpu[,pval_adj:=pval/sort(pval)[3],alpha_nuc]
  res_gpu[,method:="tshell",]
  
  return(res_gpu)
}

get_most_popularV<-function(str){
  if (length(str)!=0){
    names(sort(-table(sapply(strsplit(str,split = "*",fixed = T),"[[",1))))[1]
  }else{character()}
}

add_VJ_aa<-function(nSeqCDR3s,source_data){
  tmp<-source_data[targetSequences%in%nSeqCDR3s,.(cdr3aa=aaSeqCDR3[1],v=get_most_popularV(allVHitsWithScore),j=get_most_popularV(allJHitsWithScore)),targetSequences]
  setkey(tmp,targetSequences)
  tmp[nSeqCDR3s,]
}

add_sign<-function(tirtl_m,sem_threshold=2.5,log2FC_threshold=3,pseudo1=1e-6,pseudo2=1e-6){
  #titrl_m[,sign:=ifelse(abs(log2FC)>log2FC_threshold&sem.x<sem_threshold&sem.y<sem_threshold,"*","")]
  #column sign is "expanded" if log2FC>threshold, "contracted" is log2FC< negative threshold, "non-signif" else
  tirtl_m[,sign:="stable",]
  tirtl_m[log2FC<(-log2FC_threshold)&((avg.y+pseudo2+sem_threshold*sem.y)<(avg.x+pseudo1-sem_threshold*sem.x)),sign:="down",]
  tirtl_m[log2FC>(log2FC_threshold)&((avg.x+pseudo1+sem_threshold*sem.x)<(avg.y+pseudo2-sem_threshold*sem.y)),sign:="up",]
  return(tirtl_m)
}

merge_TIRTL<-function(mlist1,mlist2,wells1,wells2,thres1=4,thres2=4,pseudo1=1e-6,pseudo2=1e-6)
{
  tp1<-mlist1[well%in%wells1,.(avg=sum(readFraction)/length(wells1),
                               sem=sd(c(readFraction,rep(0,times=length(wells1)-length(unique(well)))))/sqrt(length(wells1)),
                               wells=length(unique(well)), 
                               avg_well=sum(readFraction)/length(unique(well))),nSeqCDR3]
  tp2<-mlist2[well%in%wells2,.(avg=sum(readFraction)/length(wells2),
                               sem=sd(c(readFraction,rep(0,times=length(wells2)-length(unique(well)))))/sqrt(length(wells2)),
                               wells=length(unique(well)),
                               avg_well=sum(readFraction)/length(unique(well))),nSeqCDR3]
  
  
  tmpm<-na_to0(merge(tp1,tp2,by="nSeqCDR3",all=T))[wells.x>thres1|wells.y>thres2,]
  
  tmpm[(wells.x<3),]$sem.x=mean(tmpm[wells.x==3,]$sem.x)*2
  tmpm[(wells.y<3),]$sem.y=mean(tmpm[wells.y==3,]$sem.y)*2
  tmpm[,log2FC:=log2((avg.y+pseudo1)/(avg.x+pseudo2)),]
  #tmpm$aaSeqCDR3=translate_cdr(tmpm$nSeqCDR3)
  #tmpm$inframe=(nchar(tmpm$nSeqCDR3)%%3)==0
  return(tmpm)
}

na_to0<-function (x) {
  x[is.na(x)]<-0
  x
}

run_longitudinal_analysis_sub<-function(folder1,folder2,well_pos1=3,well_pos2=3,well_filter_thres=0.75,wellset1=get_well_subset(1:16,1:24),wellset2=get_well_subset(1:16,1:24)) #full plates by default
{
  print("start double timepoint analysis")
  print(Sys.time())
  mlist1<-lapply(list.files(path = folder1,full.names = T),fread)
  names(mlist1)<-list.files(path = folder1,full.names = F)
  mlist1_b<-rbindlist(getb(mlist1),idcol="filename")
  mlist1_b[,well:=sapply(strsplit(filename,split="_",fixed=T),"[[",well_pos1),]
  mlist1_a<-rbindlist(geta(mlist1),idcol="filename")
  mlist1_a[,well:=sapply(strsplit(filename,split="_",fixed=T),"[[",well_pos1),]
  
  mlist2<-lapply(list.files(path = folder2,full.names = T),fread)
  names(mlist2)<-list.files(path = folder2,full.names = F)
  mlist2_b<-rbindlist(getb(mlist2),idcol="filename")
  mlist2_b[,well:=sapply(strsplit(filename,split="_",fixed=T),"[[",well_pos2),]
  mlist2_a<-rbindlist(geta(mlist2),idcol="filename")
  mlist2_a[,well:=sapply(strsplit(filename,split="_",fixed=T),"[[",well_pos2),]
  print("file read done")
  print(Sys.time())
  
  tmpm_a<-add_sign(merge_TIRTL(mlist1_a,
                               mlist2_a,
                               get_good_wells_sub(geta(mlist1),getb(mlist1),round(well_filter_thres*mean(sapply(mlist1,nrow))),pos=well_pos1,wellset = wellset1)$well_ids,
                               get_good_wells_sub(geta(mlist2),getb(mlist2),round(well_filter_thres*mean(sapply(mlist2,nrow))),pos=well_pos2,wellset = wellset2)$well_ids,4,4))
  tmpm_b<-add_sign(merge_TIRTL(mlist1_b,
                               mlist2_b,
                               get_good_wells_sub(geta(mlist1),getb(mlist1),round(well_filter_thres*mean(sapply(mlist1,nrow))),pos=well_pos1,wellset = wellset1)$well_ids,
                               get_good_wells_sub(geta(mlist2),getb(mlist2),round(well_filter_thres*mean(sapply(mlist2,nrow))),pos=well_pos2,wellset = wellset2)$well_ids,4,4))
  print("timepoint comparisons done")
  print(Sys.time())
  
  list(alpha=tmpm_a,beta=tmpm_b)
}

write_dat<-function(x,fname,rows=F){
  write.table(x,sep="\t",quote = F,row.names = rows,col.names=F,file = fname)
}

run_single_point_analysis_sub_gpu<-function(folder_path,prefix="tmp",well_filter_thres=0.5,min_reads=0,min_wells=2,well_pos=3,wellset1=get_well_subset(1:16,1:24),compute=T,backend="numpy"){ #this is with cpu backend
  print("start")
  print(Sys.time())
  mlist<-lapply(list.files(path = folder_path,full.names = T),fread)
  names(mlist)<-gsub(".","_",list.files(path = folder_path,full.names = F),fixed=T)
  mlista<-geta(mlist)
  mlistb<-getb(mlist)
  print("clonesets loaded")
  print(Sys.time())
  print(names(mlista))
  wellsub<-sapply(strsplit(names(mlista),split="_",fixed=T),"[[",well_pos)%in%wellset1
  clone_thres=round(well_filter_thres*mean(sapply(mlista,nrow)[wellsub]))
  rm(mlist)
  
  qc<-get_good_wells_sub(mlista,mlistb,clone_thres,pos=well_pos,wellset=wellset1)
  print("clone_threshold for QC:")
  print(clone_thres)
  print("alpha wells passing QC:")
  print(table(qc$a))
  print("beta wells passing QC:")
  print(table(qc$b))
  
  mlista<-mlista[qc$a]#downsize to qc
  mlistb<-mlistb[qc$b]#downsize to qc
  
  print(Sys.time())
  print("Merging alpha clonesets...")
  bigma<-big_merge_freqs2(mlista,min_reads = min_reads)
  print("Done! Unique alpha clones and wells after filtering:")
  print(dim(bigma))
  bigmas<-bigma[rowSums(bigma>0)>min_wells,]
  print(Sys.time())
  print(paste0("Unique alpha clones and wells in more than: ",min_wells," wells"))
  print(dim(bigmas))
  print("Merging beta clonesets...")
  bigmb<-big_merge_freqs2(mlistb,min_reads = min_reads)
  print("Done! beta clones and wells after filtering:")
  print(dim(bigmb))
  bigmbs<-bigmb[rowSums(bigmb>0)>min_wells,]
  print(Sys.time())
  print(paste0("Unique beta clones and wells in more than: ",min_wells," wells"))
  print(dim(bigmbs))
  print("Writing files for back-end pairing script...")
  write(rownames(bigmas),file=paste0(prefix,"_bigmas_names.tsv"))
  write(rownames(bigmbs),file=paste0(prefix,"_bigmbs_names.tsv"))
  write_dat(as.matrix(bigmas),fname = paste0(prefix,"_bigmas.tsv"))
  write_dat(as.matrix(bigmbs),fname = paste0(prefix,"_bigmbs.tsv"))
  print(Sys.time())
  n_wells=ncol(bigmas)
  print("Pre-computing look-up table:")
  mdh<-madhyper_surface(n_wells = ncol(bigmas),cells = clone_thres,alpha=2,prior = 1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5)
  write_dat(mdh,fname = paste0(prefix,"_mdh.tsv"))
  print(Sys.time())
  if(compute==T)
    if(backend=="cupy")system(paste0("python3 cupy_backend_script.py ",prefix,collapse=""))
    else if(backend=="mlx")system(paste0("python3 mlx_backend_script.py ",prefix,collapse=""))
      else{system(paste0("python3 numpy_backend_script.py ",prefix,collapse=""))}

  print("Loading and filtering results, adding amino acid and V segment information")
  
    # and here we go read it: 
  gpu_res<-read_gpu(prefix)
  gpu_res_corr<-read_gpu_corr(prefix)
  # I also want to compute 
  result<-rbind(gpu_res,gpu_res_corr,fill=T)

  result[,loss_a_frac:=(wb-wij)/(wij+(wb-wij)+(wa-wij)),]
  result[,loss_b_frac:=(wa-wij)/(wij+(wb-wij)+(wa-wij)),]
  result[,wi:=wa-wij,]
  result[,wj:=wb-wij,]
  
  #groom the output. 
  
  unique_combinations <- unique(result[, .(wi, wj, wij)])
  for (i in 1:nrow(unique_combinations)){
    #if(i%%1000==0)print(i)
    unique_combinations$score[i]<-log10(estimate_pair_prob(wi = unique_combinations[i,]$wi,wj = unique_combinations[i,]$wj,w_ij = unique_combinations[i,]$wij,n_wells,cpw = clone_thres,alpha = 2,prior=1/(as.numeric(nrow(bigmas))*(as.numeric(nrow(bigmbs))))**0.5))
  }
  
  #result<-result[order(-method),][!duplicated(alpha_beta),]#version without filter
  result<-merge(result, unique_combinations, by = c("wi", "wj", "wij"), all.x = TRUE)[method=="madhype"|(method=="tshell"&wij>2&pval_adj<1e-10&(loss_a_frac+loss_b_frac)<0.5),]
  
  #result<-result[(((loss_a_frac+loss_b_frac)<0.5)&(wij>3))|(score>0.1),]
  
  tp_a<-add_VJ_aa(result$alpha_nuc,rbindlist(mlista))
  result$cdr3a=tp_a$cdr3aa
  result$va=tp_a$v
  result$ja=tp_a$j
  
  tp_b<-add_VJ_aa(result$beta_nuc,rbindlist(mlistb))
  result$cdr3b=tp_b$cdr3aa
  result$vb=tp_b$v
  result$jb=tp_b$j
  print(Sys.time())
  print("All is done! Number of paired clones:")
  print(table(result$method))
  fwrite(result,paste0(prefix,"_TIRTLoutput.tsv"),sep="\t")    
  return(result)
}

get_clonotypes_10x<-function(TCRs){ #this makes neat table from filtered contig annotations.
  ctg<-data.table(V1=unique(TCRs$barcode))
  ctg$cdr3b<-NA
  ctg$cdr3b_nt<-NA
  ctg$vb<-NA
  ctg$jb<-NA
  ctg$cdr3a<-NA
  ctg$cdr3a_nt<-NA
  ctg$va<-NA
  ctg$ja<-NA
  ctg$cdr3a2<-NA
  ctg$cdr3a2_nt<-NA
  ctg$va2<-NA
  ctg$ja2<-NA
  for (i in 1:nrow(ctg)){
    if (i%%1000==0) print(i)
    ctg$cdr3b[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$cdr3
    ctg$cdr3b_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$cdr3_nt
    ctg$vb[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$v_gene
    ctg$jb[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRB",,][1,,]$j_gene
    ctg$cdr3a[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$cdr3
    ctg$cdr3a_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$cdr3_nt
    ctg$va[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$v_gene
    ctg$ja[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][1,,]$j_gene
    ctg$cdr3a2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$cdr3
    ctg$cdr3a2_nt[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$cdr3_nt
    ctg$va2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$v_gene
    ctg$ja2[i]<-TCRs[barcode==ctg$V1[i],,][order(-umis)][chain=="TRA",,][2,,]$j_gene
  }
  ctg[,n_cells:=.N,.(cdr3b_nt,cdr3a_nt)]#N cells/clone
  ctg
}


process_10x<-function(path){
  dt_10x<-fread(path)
  dt_10x_clean<-get_clonotypes_10x(dt_10x)
  dt_10x_complete<-dt_10x_clean[order(-n_cells),][!duplicated(paste0(cdr3b_nt,"_",cdr3a_nt)),][!is.na(cdr3b_nt)&!is.na(cdr3a_nt),] 
  dt_10x_complete[,alpha_beta:=paste0(cdr3a_nt,"_",cdr3b_nt),]
  dt_10x_complete[,beta_nuc:=cdr3b_nt,]
  dt_10x_complete[,alpha_nuc:=cdr3a_nt,]
  list(complete=dt_10x_complete,clean=dt_10x_clean,raw=dt_10x)
}


#test<-run_single_point_analysis_sub_gpu("preprint_data_combined/MVP093_PAIRSEQ/plate1/",well_pos = 4,wellset1 = get_well_subset(1:16,19:24))






