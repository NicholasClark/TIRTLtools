#' Re-order samples in a TIRTLseqDataSet object
#'
#' `r lifecycle::badge('experimental')`
#' This function changes the order of the samples in a TIRTLseqDataSet object to the order that the user specifies.
#'
#' @param data a TIRTLseqDataSet object
#' @param samples the samples of 'data', in their desired order.
#'
#' Either:
#' 1. A numeric vector of indices
#' 2. A character vector of sample names
#' 3. A character vector of metadata conditions where each entry is of the form \code{"column==value"}
#'
#' @return A TIRTLseqDataSet object with re-ordered samples.
#'
#' @family data_wrangling
#'
#' @examples
#' load_example_data(dataset = "SJTRC_longitudinal")
#' print(SJTRC_longitudinal$meta)
#' new_order = names(SJTRC_longitudinal$data) %>% rev()
#' print(new_order)
#' SJTRC_longitudinal = reorder_samples(SJTRC_longitudinal, new_order)
#' print(SJTRC_longitudinal$meta)
#'
reorder_samples = function(data, samples) {
  if(length(samples) != length(data$data)) stop("'samples' must be the same length as the number of samples in 'data'.")
  filter_dataset(data, samples)
}
