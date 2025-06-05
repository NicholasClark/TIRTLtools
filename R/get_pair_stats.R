get_pair_stats = function(data) {
  pair_stats = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    sample_name = names(data$data)[i]
    print(i)
    df = .get_pair_stats_single(x) %>% mutate(sample_id = sample_name)
  }) %>% bind_rows()
  return(pair_stats)
}
