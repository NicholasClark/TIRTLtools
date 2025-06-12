reorder_samples = function(data, samples) {
  if(length(samples) != length(data$data)) stop()
  filter_dataset(data, samples)
}
