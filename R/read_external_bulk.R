#' Read and process bulk single-chain TCR-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' 
#' This function reads and processes bulk single-chain TCR-sequencing data from a non-TIRTLseq assay.
#' Currently, only MiXCR format is supported.
#'
#' @details
#' Supported data types:
#' - \code{"MiXCR"} - "<sample_name>_TRA.tsv" and "<sample_name>_TRB.tsv"
#'
#' @param path the path to the data file
#' @param format the data format. Currently only "MiXCR" is supported.
#'
#' @return A list containing the following slots:
#' - df - a data frame with a few columns modified and renamed
#' - df_raw - the original data, un-modified
#'
#' @family data_processing

read_external_bulk = function(path, format = "MiXCR") {
  df_orig = fread(path)
  cols = .get_necessary_cols_bulk()
  .check_cols(df_orig, cols) ## errors if it doesn't have the needed columns
  df = df_orig %>%
    mutate(v = gsub("\\*.*", "", allVHitsWithScore),
           j = gsub("\\*.*", "", allJHitsWithScore)) %>%
    select(targetSequences, readCount, v, j, aaSeqCDR3, readFraction, everything())
  return(list(df = df, df_raw = df_orig))
}
