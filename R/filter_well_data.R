## input is a data frame with data from all wells
filter_well_data = function(well_data_df, min_fraction = 1e-6) {
  well_data_df = well_data_df %>%
    mutate(v = get_segments_all(allVHitsWithScore), j = get_segments_all(allJHitsWithScore))  %>%
    mutate(id = paste(targetSequences, v, j, aaSeqCDR3))
  pseudobulk = well_data_df %>%
    group_by(targetSequences, v, j, aaSeqCDR3, id) %>%
    summarize(readCount_total = sum(readCount),
              allVHitsWithScore = allVHitsWithScore[1],
              allJHitsWithScore = allJHitsWithScore[1]
    ) %>%
    ungroup() %>%
    mutate(vj = paste(v,j)) %>%
    mutate(readFraction_total = readCount_total/sum(readCount_total)) %>%
    arrange(desc(readFraction_total))

  pseudobulk_sub = pseudobulk %>% filter(readFraction_total >= min_fraction)

  well_data_df_sub = well_data_df %>% filter(id %in% pseudobulk_sub$id) %>%
    select(-id)
  return(well_data_df_sub)
}
