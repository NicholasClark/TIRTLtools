
#TCRs = fread("~/git/newell_benchmarking/data/GSM8951948_GEMX5P_F5_filtered_contig_annotations.csv")

process_10x_old<-function(path){
  dt_10x<-fread(path)
  dt_10x_clean<-get_clonotypes_10x_old(dt_10x)
  dt_10x_complete<-dt_10x_clean[order(-n_cells),][!duplicated(paste0(cdr3b_nt,"_",cdr3a_nt)),][!is.na(cdr3b_nt)&!is.na(cdr3a_nt),]
  dt_10x_complete[,alpha_beta:=paste0(cdr3a_nt,"_",cdr3b_nt),]
  dt_10x_complete[,beta_nuc:=cdr3b_nt,]
  dt_10x_complete[,alpha_nuc:=cdr3a_nt,]
  list(complete=dt_10x_complete,clean=dt_10x_clean,raw=dt_10x)
}


get_clonotypes_10x_old<-function(TCRs){ #this makes neat table from filtered contig annotations.
  ## make data frame with one row for each unique barcode
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
  #' steps:
  #' 1) Filter data to each unique barcode, for one chain (alpha or beta)
  #' 2) Order by most to least UMIs
  #' 3) Record v,j, cdr3_aa, cdr3_nt for top beta and top two alphas.
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
