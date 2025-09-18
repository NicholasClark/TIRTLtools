#' Change the order of samples in a TIRTLseqData object
#'
#' @family data_wrangling
#'
reorder_samples = function(data, samples) {
  if(length(samples) != length(data$data)) stop()
  filter_dataset(data, samples)
}
