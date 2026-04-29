#' Process single-cell TIRTL-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function reads a directory of MiXCR output files from single-cell TIRTL-seq
#' and outputs a data frame with a single TCR for each well.
#'
#' @param folder the path of the folder with .tsv (or .tsv.gz) output files from MiXCR
#' @param min_read_fraction the minimum read fraction required to select a chain for a well
#' @param min_read_count the minimum read count required to select a chain for a well
#' @param require_functional whether to require a functional chain (i.e. without frameshift or stop codons in the CDR3 region)
#'
#' @return
#' A data frame with one TCR pair selected per well.
#' @family single-cell
#' @export

process_scTIRTLseq = function(folder, min_read_fraction = 0.1, min_read_count = 15, require_functional = TRUE){
#read_all_files_TIRTLseq<-function(folder){
  file_list <- list.files(path = folder, pattern = "(TCRa|TCRb).*\\.(tsv|tsv\\.gz)$", full.names = TRUE)
  #file_list<-gsub(".","_",file_list,fixed=T)
  #file_list<-gsub("TRA","TRAD",file_list)
  # Read and process all the files
  #print(file_list)
  message(paste("Loading", length(file_list), "files..."))
  processed_files <- lapply(1:length(file_list), function(i) {
    x=file_list[i]
    tmp = read_and_process_file_TIRTLseq(x,
                                   min_read_fraction = min_read_fraction,
                                   min_read_count = min_read_count,
                                   require_functional = require_functional
                                   )
    if( i %% 100 == 0 | i %% length(file_list) == 0 ) message(paste(i, "files loaded"))
    return(tmp)
  })
  # Combine the results into a single table
  combined_dt <- rbindlist(processed_files)
  aggr_fun <- function(x) {
    if (length(x) == 0) {
      return(NA)
    } else if (length(x) == 1) {
      return(x)
    } else {
      return(x[1])
    }
  }

  # Spread the results into separate columns for TRB and TRAD
  result <- data.table::dcast(combined_dt, barcode + x + y ~ chain_type + alpha_no, value.var = c("nSeqCDR3", "aaSeqCDR3", "allVHitsWithScore", "allJHitsWithScore", "cloneId"), fun.aggregate = aggr_fun)

  # the following is transforming to the legacy file format
  data.table::setnames(result,
           old = c("nSeqCDR3_clones_TRA_first", "nSeqCDR3_clones_TRA_second", "nSeqCDR3_clones_TRB_first", "aaSeqCDR3_clones_TRA_first", "aaSeqCDR3_clones_TRA_second", "aaSeqCDR3_clones_TRB_first", "allVHitsWithScore_clones_TRA_first", "allVHitsWithScore_clones_TRA_second", "allVHitsWithScore_clones_TRB_first", "allJHitsWithScore_clones_TRA_first", "allJHitsWithScore_clones_TRA_second", "allJHitsWithScore_clones_TRB_first", "cloneId_clones_TRA_first", "cloneId_clones_TRA_second", "cloneId_clones_TRB_first"),
           new = c("nSeqCDR3_alpha", "nSeqCDR3_alpha_second", "nSeqCDR3_beta", "aaSeqCDR3u_alpha", "aaSeqCDR3u_alpha_second", "aaSeqCDR3u_beta", "bestV_alpha", "bestV_alpha_second", "bestV_beta", "bestJ_alpha", "bestJ_alpha_second", "bestJ_beta", "clone_rank_alpha", "clone_rank_alpha_second", "clone_rank_beta"))

  reordered_result <- result[, .(plate = NA_character_,
                                 barcode,
                                 read = NA_character_,  # Add an empty column for 'read_beta'
                                 nSeqCDR3_beta,
                                 aaSeqCDR3u_beta,
                                 bestV_beta=splitvj(bestV_beta),
                                 bestJ_beta=splitvj(bestJ_beta),
                                 TCR_rank_beta = clone_rank_beta,
                                 N_beta = NA_integer_,  # Add an empty column for 'N_beta'
                                 N_prop_beta = NA_real_,  # Add an empty column for 'N_prop_beta'
                                 clone_rank_beta,
                                 functional_beta = NA,  # Add an empty column for 'functional_beta'
                                 read_alpha = NA_character_,  # Add an empty column for 'read_alpha'
                                 nSeqCDR3_alpha,
                                 aaSeqCDR3u_alpha,
                                 bestV_alpha=splitvj(bestV_alpha),
                                 bestJ_alpha=splitvj(bestJ_alpha),
                                 TCR_rank_alpha = clone_rank_alpha,
                                 N_alpha = NA_integer_,  # Add an empty column for 'N_alpha'
                                 N_prop_alpha = NA_real_,  # Add an empty column for 'N_prop_alpha'
                                 clone_rank_alpha,
                                 functional_alpha = NA,  # Add an empty column for 'functional_alpha'
                                 read_alpha_second = NA_character_,  # Add an empty column for 'read_alpha_second'
                                 nSeqCDR3_alpha_second,
                                 aaSeqCDR3u_alpha_second,
                                 bestV_alpha_second=splitvj(bestV_alpha_second),
                                 bestJ_alpha_second=splitvj(bestJ_alpha_second),
                                 TCR_rank_alpha_second = clone_rank_alpha_second,
                                 N_alpha_second = NA_integer_,  # Add an empty column for 'N_alpha_second'
                                 N_prop_alpha_second = NA_real_,  # Add an empty column for 'N_prop_alpha_second'
                                 clone_rank_alpha_second,
                                 functional_alpha_second = NA,  # Add an empty column for 'functional_alpha_second'
                                 bestN_alpha_second = NA_character_,  # Add an empty column for 'bestN_alpha_second'
                                 x,
                                 y)]
  reordered_result
}


mark_functional<-function(dt){
  dt[,functional:=T,]
  dt[aaSeqCDR3=="",functional:=F,]
  dt[grepl("_",aaSeqCDR3,fixed=T),functional:=F,]
  dt[grepl("*",aaSeqCDR3,fixed=T),functional:=F,]
  dt
}


read_and_process_file_TIRTLseq<-function(file_path, min_read_fraction = 0.1, min_read_count = 15, require_functional = TRUE) {
  dt <- fread(file_path)
  # Extract the well information and chain type from the file name
  #well_info <- unlist(strsplit(basename(file_path), "\\."))[2]
  well_info = .get_well_from_filename(file_path)
  #print(well_info)
  chain_type <- unlist(strsplit(basename(file_path), "\\."))[3]
  dt[, c("barcode", "x", "y", "chain_type") := .(well_info,  substr(well_info, 1, 1), as.integer(substr(well_info, 2, nchar(well_info))), chain_type)]
  dt <- mark_functional(dt)
  dt <- dt[(readFraction >= min_read_fraction) & (readCount >= min_read_count) & (functional == require_functional)] # filter here! # was >50
  # If more than one row and chain_type is "TRA", sort by readFraction and take the first two rows
  if (chain_type == "clones_TRA" && nrow(dt) > 1) {
    dt <- dt[order(-readFraction)][1:2, ]
    dt[, alpha_no := c("first", "second")]
  } else {
    dt[, alpha_no := "first"]
  }
  return(dt)
}


splitvj<-function(x){
  sapply(strsplit(x,"*",fixed=T),"[[",1)
}
