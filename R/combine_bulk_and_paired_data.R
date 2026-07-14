#' Combine bulk single-chain TCR data with paired TCR data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' 
#' This function...
#'
#' @param alpha a dataframe with alpha chain bulk TCR data
#' @param beta a dataframe with beta chain bulk TCR data
#' @param paired a dataframe with paired chain TCR data
#' @param meta either the sample name or a metadata data frame containing the sample name in the 'sample_id' column.
#'
#' @return a pairedTCRDataSet object
#'
#' @family data_processing
#'
combine_bulk_and_paired_data = function(alpha, beta, paired, meta) {
  assert_data_frame(alpha)
  assert_data_frame(beta)
  assert_data_frame(paired)
  assert_subset(c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "readFraction"), colnames(alpha))
  assert_subset(c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "readFraction"), colnames(beta))
  assert_subset(c("va", "ja", "alpha_nuc", "beta_nuc", "cdr3a", "cdr3b"), colnames(paired))
  is_df_meta = checkmate::test_data_frame(meta)
  is_char_meta = test_character(meta)
  if(is_df_meta) {
    assert(nrow(meta) == 1)
    assert_subset("sample_id", colnames(meta))
    sample_name = meta$sample_id
  } else if(is_char_meta) {
    assert(length(meta) == 1)
    sample_name = meta
    meta = tibble(sample_id = sample_name)
  } else {
    stop("'meta' must be either the sample name or a metadata data frame")
  }
  obj = list(alpha = alpha, beta = beta, paired = paired)
  obj = obj %>% .process_TIRTLseq_single(clean = FALSE, remove_nonfunctional = FALSE, verbose = FALSE)
  class(obj) = "pairedTCRData"
  data_list = list(obj)
  names(data_list) = sample_name
  out = list(meta = meta, data = data_list)
  class(out) = "pairedTCRDataSet"
  return(out)
}
