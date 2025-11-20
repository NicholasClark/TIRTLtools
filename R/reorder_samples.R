#' Re-order samples in a TIRTLseqData object
#'
#' `r lifecycle::badge('experimental')`
#' This function changes the order of the samples in a TIRTLseqData object to the order that the user specifies.
#'
#' @param data a TIRTLseqData object
#' @param samples the samples of 'data', in their desired order.
#'
#' Either:
#' 1. A numeric vector of indices
#' 2. A character vector of sample names
#' 3. A character vector of metadata conditions where each entry is of the form \code{"column==value"}
#'
#' @return A TIRTLseqData object with re-ordered samples.
#'
#' @family data_wrangling
#'
reorder_samples = function(data, samples) {
  if(length(samples) != length(data$data)) stop("'samples' must be the same length as the number of samples in 'data'.")
  filter_dataset(data, samples)
}
