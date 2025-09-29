#' Change the order of samples in a TIRTLseqData object
#'
#' @param data a TIRTLseqData object
#' @param samples the samples of 'data', in their desired order. Either 1) a numeric vector of indices,
#' 2) a character vector of sample names, or 3) a character vector of metadata
#' conditions where each entry is of the form \code{"column==value"}.
#'
#' @family data_wrangling
#'
reorder_samples = function(data, samples) {
  if(length(samples) != length(data$data)) stop("'samples' must be the same length as the number of samples in 'data'.")
  filter_dataset(data, samples)
}
