# library(data.table)
# read_all_files_TIRTLseq<-function(folder){
#   file_list <- list.files(path = folder, pattern = "(TCRa|TCRb).*\\.(tsv|tsv\\.gz)$", full.names = TRUE)
#   #file_list<-gsub(".","_",file_list,fixed=T)
#   #file_list<-gsub("TRA","TRAD",file_list)
#   # Read and process all the files
#   #print(file_list)
#   processed_files <- lapply(file_list, read_and_process_file_TIRTLseq)
#   # Combine the results into a single table
#   combined_dt <- rbindlist(processed_files)
#   aggr_fun <- function(x) {
#     if (length(x) == 0) {
#       return(NA)
#     } else if (length(x) == 1) {
#       return(x)
#     } else {
#       return(x[1])
#     }
#   }
#
#   # Spread the results into separate columns for TRB and TRAD
#   result <- dcast(combined_dt, barcode + x + y ~ chain_type + alpha_no, value.var = c("nSeqCDR3", "aaSeqCDR3", "allVHitsWithScore", "allJHitsWithScore", "cloneId"), fun.aggregate = aggr_fun)
#
#   # the following is transforming to the legacy file format
#   setnames(result,
#            old = c("nSeqCDR3_clones_TRA_first", "nSeqCDR3_clones_TRA_second", "nSeqCDR3_clones_TRB_first", "aaSeqCDR3_clones_TRA_first", "aaSeqCDR3_clones_TRA_second", "aaSeqCDR3_clones_TRB_first", "allVHitsWithScore_clones_TRA_first", "allVHitsWithScore_clones_TRA_second", "allVHitsWithScore_clones_TRB_first", "allJHitsWithScore_clones_TRA_first", "allJHitsWithScore_clones_TRA_second", "allJHitsWithScore_clones_TRB_first", "cloneId_clones_TRA_first", "cloneId_clones_TRA_second", "cloneId_clones_TRB_first"),
#            new = c("nSeqCDR3_alpha", "nSeqCDR3_alpha_second", "nSeqCDR3_beta", "aaSeqCDR3u_alpha", "aaSeqCDR3u_alpha_second", "aaSeqCDR3u_beta", "bestV_alpha", "bestV_alpha_second", "bestV_beta", "bestJ_alpha", "bestJ_alpha_second", "bestJ_beta", "clone_rank_alpha", "clone_rank_alpha_second", "clone_rank_beta"))
#
#   reordered_result <- result[, .(plate = NA_character_,
#                                  barcode,
#                                  read = NA_character_,  # Add an empty column for 'read_beta'
#                                  nSeqCDR3_beta,
#                                  aaSeqCDR3u_beta,
#                                  bestV_beta=splitvj(bestV_beta),
#                                  bestJ_beta=splitvj(bestJ_beta),
#                                  TCR_rank_beta = clone_rank_beta,
#                                  N_beta = NA_integer_,  # Add an empty column for 'N_beta'
#                                  N_prop_beta = NA_real_,  # Add an empty column for 'N_prop_beta'
#                                  clone_rank_beta,
#                                  functional_beta = NA,  # Add an empty column for 'functional_beta'
#                                  read_alpha = NA_character_,  # Add an empty column for 'read_alpha'
#                                  nSeqCDR3_alpha,
#                                  aaSeqCDR3u_alpha,
#                                  bestV_alpha=splitvj(bestV_alpha),
#                                  bestJ_alpha=splitvj(bestJ_alpha),
#                                  TCR_rank_alpha = clone_rank_alpha,
#                                  N_alpha = NA_integer_,  # Add an empty column for 'N_alpha'
#                                  N_prop_alpha = NA_real_,  # Add an empty column for 'N_prop_alpha'
#                                  clone_rank_alpha,
#                                  functional_alpha = NA,  # Add an empty column for 'functional_alpha'
#                                  read_alpha_second = NA_character_,  # Add an empty column for 'read_alpha_second'
#                                  nSeqCDR3_alpha_second,
#                                  aaSeqCDR3u_alpha_second,
#                                  bestV_alpha_second=splitvj(bestV_alpha_second),
#                                  bestJ_alpha_second=splitvj(bestJ_alpha_second),
#                                  TCR_rank_alpha_second = clone_rank_alpha_second,
#                                  N_alpha_second = NA_integer_,  # Add an empty column for 'N_alpha_second'
#                                  N_prop_alpha_second = NA_real_,  # Add an empty column for 'N_prop_alpha_second'
#                                  clone_rank_alpha_second,
#                                  functional_alpha_second = NA,  # Add an empty column for 'functional_alpha_second'
#                                  bestN_alpha_second = NA_character_,  # Add an empty column for 'bestN_alpha_second'
#                                  x,
#                                  y)]
#   reordered_result
# }
#
#
# mark_functional<-function(dt){
#   dt[,functional:=T,]
#   dt[aaSeqCDR3=="",functional:=F,]
#   dt[grepl("_",aaSeqCDR3,fixed=T),functional:=F,]
#   dt[grepl("*",aaSeqCDR3,fixed=T),functional:=F,]
#   dt
# }
#
#
# read_and_process_file_TIRTLseq<-function(file_path) {
#   dt <- fread(file_path)
#   # Extract the well information and chain type from the file name
#   #well_info <- unlist(strsplit(basename(file_path), "\\."))[2]
#   well_info = .get_well_from_filename(file_path)
#   #print(well_info)
#   chain_type <- unlist(strsplit(basename(file_path), "\\."))[3]
#   dt[, c("barcode", "x", "y", "chain_type") := .(well_info,  substr(well_info, 1, 1), as.integer(substr(well_info, 2, nchar(well_info))), chain_type)]
#   dt <- mark_functional(dt)
#   dt <- dt[(readFraction >= 0.1) & (readCount>15) & (functional == TRUE)] # filter here! # was >50
#   # If more than one row and chain_type is "TRA", sort by readFraction and take the first two rows
#   if (chain_type == "clones_TRA" && nrow(dt) > 1) {
#     dt <- dt[order(-readFraction)][1:2, ]
#     dt[, alpha_no := c("first", "second")]
#   } else {
#     dt[, alpha_no := "first"]
#   }
#   return(dt)
# }
#
#
# splitvj<-function(x){
#   sapply(strsplit(x,"*",fixed=T),"[[",1)
# }
#
#
#
# library(ggplot2)
# library(ggplate)
# library(patchwork)
# TIRTL_pallet2<-c("#007EA7","#3CAF82","#FFD16E","#A53828","#F37748","grey90")
# substitution3 <- c("TRUETRUE" = "ab", "TRUEFALSE" = "a lost", "FALSETRUE" = "b lost")
#
#
# #Tara_plate_MVP158_UD05<-read_all_files_TIRTLseq("/Volumes/TIRTL/tirtlSeq_data/MVP158_DNMT3_RTstorage_TaraTIL/clones_MVP158_DNMT3_RTstorage/TCR_clones_UD05/")
# #Tara_plate_MVP158_UD06<-read_all_files_TIRTLseq("/Volumes/TIRTL/tirtlSeq_data/MVP158_DNMT3_RTstorage_TaraTIL/clones_MVP158_DNMT3_RTstorage/TCR_clones_UD06/")
# noroplates<-list(NoroVirus_Exp_NIH3_Healthy_single=read_all_files_TIRTLseq("~/git/from_Misha/clones_MVP224_Norovirus_expansion/TCR_clones_UD17"),#, well_filter_thres=0.5,wellset1 = get_well_subset(1:16,1:24),pseudobulk=T,backend="cupy",prefix="NoroVirus_Exp_NIH3_Healthy_single")
#                  NoroVirus_Exp_NIH1_single=read_all_files_TIRTLseq("~/git/from_Misha/clones_MVP224_Norovirus_expansion/TCR_clones_UD19"),#, well_filter_thres=0.5,wellset1 = get_well_subset(1:16,1:24),pseudobulk=T,backend="cupy",prefix="NoroVirus_Exp_NIH1_single")
#                  NoroVirus_Exp_NIH2_single=read_all_files_TIRTLseq("~/git/from_Misha/clones_MVP224_Norovirus_expansion/TCR_clones_UD21"),#, well_filter_thres=0.5,wellset1 = get_well_subset(1:16,1:24),pseudobulk=T,backend="cupy",prefix="NoroVirus_Exp_NIH2_single")
#                  NoroVirus_Exp_NIH4_single=read_all_files_TIRTLseq("~/git/from_Misha/clones_MVP224_Norovirus_expansion/TCR_clones_UD23"),#, well_filter_thres=0.5,wellset1 = get_well_subset(1:16,1:24),pseudobulk=T,backend="cupy",prefix="NoroVirus_Exp_NIH4_single")
#                  NoroVirus_Exp_NIH5_single=read_all_files_TIRTLseq("~/git/from_Misha/clones_MVP224_Norovirus_expansion/TCR_clones_UD25"))#, well_filter_thres=0.5,wellset1 = get_well_subset(1:16,1:24),pseudobulk=T,backend="cupy",prefix="NoroVirus_Exp_NIH5_250")
# noroplatesb<-rbindlist(noroplates,idcol = "plate")
# fwrite(noroplatesb,sep="\t",file="~/git/from_Misha/norovirus_one_cell_well.tsv")
#
# i=3
# plate_plot(noroplates[[i]][,.(barcode,chains=substitution3[paste0(!is.na(nSeqCDR3_beta),!is.na(nSeqCDR3_alpha))]),],
#            plate_size = 384,
#            position = barcode,
#            value=chains,show_legend = FALSE,title=NULL)+scale_fill_manual(values=c("ab"=TIRTL_pallet2[2],"b lost"=TIRTL_pallet2[5],"a lost"=TIRTL_pallet2[3]))+ggtitle(names(noroplates)[[i]])
#
#
# ### test with data from MVP303
#
# dir1 = "/Volumes/thomas_p/grp/TIRTL/projects/MVP303/clones/UD10"
# test = read_all_files_TIRTLseq(dir1)
#
# plate_plot(test[,.(barcode,chains=substitution3[paste0(!is.na(nSeqCDR3_beta),!is.na(nSeqCDR3_alpha))]),],
#            plate_size = 384,
#            position = barcode,
#            value=chains,show_legend = FALSE,title=NULL)+scale_fill_manual(values=c("ab"=TIRTL_pallet2[2],"b lost"=TIRTL_pallet2[5],"a lost"=TIRTL_pallet2[3]))+ggtitle("testing")
# #
# # #twoplates<-rbind(aph56_plate10,aph56_plate11)
# # completenoro<-noroplatesb[plate=="NoroVirus_Exp_NIH2_single"&!is.na(nSeqCDR3_beta)&!is.na(nSeqCDR3_alpha),]
# # completenoro[,alpha_beta:=paste0(nSeqCDR3_alpha,"_",nSeqCDR3_beta),]
# # completenoro[,alpha_beta2:=paste0(nSeqCDR3_alpha_second,"_",nSeqCDR3_beta),]
# # completenoro[,n_cells:=.N,alpha_beta]
# # completenoro[,alpha_nuc:=nSeqCDR3_alpha,]
# # completenoro[,beta_nuc:=nSeqCDR3_beta,]
# # completenoro<-completenoro[!duplicated(alpha_beta),]
# # nih2<-fread("~/git/norovirus_tirtlseq/NoroVirus_Exp_NIH2_250_TIRTLoutput.tsv")
# # concordance(completenoro,nih2)
# # concordance(complete,exp1_6plate)
# # concordance(complete,allplates)
# #
# # p15<-plate_plot(Tara_plate_MVP158_UD05[,.(barcode,chains=substitution3[paste0(!is.na(nSeqCDR3_beta),!is.na(nSeqCDR3_alpha))]),],
# #                 plate_size = 384,
# #                 position = barcode,
# #                 value=chains,show_legend = FALSE,title=NULL)+scale_fill_manual(values=c("ab"=TIRTL_pallet2[2],"b lost"=TIRTL_pallet2[5],"a lost"=TIRTL_pallet2[3]))
# # p16<-plate_plot(Tara_plate_MVP158_UD06[,.(barcode,chains=substitution3[paste0(!is.na(nSeqCDR3_beta),!is.na(nSeqCDR3_alpha))]),],
# #                 plate_size = 384,
# #                 position = barcode,
# #                 value=chains,show_legend = FALSE,title=NULL)+scale_fill_manual(values=c("ab"=TIRTL_pallet2[2],"b lost"=TIRTL_pallet2[5],"a lost"=TIRTL_pallet2[3]))
# #
# # p15+p16
# #ggsave("plate15_16.pdf",width=15,height=8)
# #Tara_plate_MVP158_UD05[barcode%in%c("C23","P19","J24","N6","N11","P8"),DNMT_mut:=TRUE,]
# #Tara_plate_MVP158_UD06[barcode%in%c("H6","I22","M24","G13"),DNMT_mut:=TRUE,]
#
# # twoplates<-rbindlist(list(plate_MVP158_UD05=Tara_plate_MVP158_UD05,plate_MVP158_UD06=Tara_plate_MVP158_UD06),idcol="plate")
# #
# # fwrite(twoplates,sep="\t",file = "Tara_plate_MVP158.tsv")
