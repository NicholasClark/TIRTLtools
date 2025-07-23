
#' input: well_data is a list of data frames with well data, each with columns:
#' "readCount", "targetSequences", "well", "allVHitsWithScore", "allJHitsWithScore"
#' function only consolidates clones with overall readFraction >= min_fraction
#' output: a data frame with data for all wells (labeled by well), where similar clones have
#' been aggregated.
consolidate_clones = function(well_data, min_fraction = -1, dist_method = "lv", linkage = "single", cut = 2) {
  well_data = bind_rows(well_data)
  well_data = well_data %>%
    mutate(v = get_all_vs(allVHitsWithScore), j = get_all_vs(allJHitsWithScore))
  pseudobulk = well_data %>%
    group_by(targetSequences, v, j) %>%
    summarize(readCount_total = sum(readCount)) %>%
    ungroup() %>%
    mutate(vj = paste(v,j)) %>%
    mutate(readFraction_total = readCount_total/sum(readCount_total)) %>%
    arrange(desc(readFraction_total))

  pseudobulk_sub = pseudobulk %>% filter(readFraction_total >= min_fraction)
  pseudobulk_extra = pseudobulk %>% filter(readFraction_total <= min_fraction) %>%
    mutate(cluster_name = paste("not_consolidated", seq_along(readFraction_total), sep = "_"), consensusClone = targetSequences)
  tic()
  pseudobulk_grp = pseudobulk_sub %>% group_by(vj) %>%
    mutate(cluster_name = paste(v,j,
      cluster_clones_hclust(targetSequences, dist_method = dist_method, linkage = linkage, cut = cut), sep = "_")) %>%
    ungroup()
  cluster_rep = pseudobulk_grp %>%
    group_by(cluster_name) %>%
    slice_max(order_by = readCount_total, n = 1, with_ties = FALSE)
  cluster_rep = cluster_rep %>% select(cluster_name, targetSequences) %>%
    dplyr::rename(consensusClone = targetSequences)
  pseudobulk_grp = pseudobulk_grp %>% left_join(cluster_rep, by = "cluster_name")
  pseudobulk_grp = bind_rows(pseudobulk_grp, pseudobulk_extra)
  toc()
  join_df = pseudobulk_grp %>% select(v,j,targetSequences, cluster_name, consensusClone)
  well_data = well_data %>% left_join(join_df, by = c("v","j", "targetSequences"))
  well_data_fold = well_data %>% group_by(consensusClone, cluster_name, v, j, well, chain, well_chain) %>%
    summarize(
      all_clones_nt = paste(targetSequences, collapse = "|"),
      all_readCount = paste(readCount, collapse = "|"),
      all_readFraction = paste(readFraction, collapse = "|"),
      readCount_max = max(readCount),
      readFraction_max = max(readFraction),
      readCount = sum(readCount),
      readFraction = sum(readFraction)) %>%
    ungroup() %>%
    arrange(desc(readCount))

  return(well_data_fold)
}

# input is a vector of nucleotide sequences for TCRs
# output is a vector of group numbers, all in the same group have Levenshtein distance <= 2 to some clone
cluster_clones_hclust = function(clones, dist_method = "lv", linkage = "single", cut = 2) {
  if(length(clones) == 1) return(1L)
  d = stringdistmatrix(clones, clones, method = dist_method)
  hc = hclust(as.dist(d), method = linkage)
  grp = cutree(hc, h = cut)
  return(grp)
}

# input is a string from the "AllVHitsWithScore" column, e.g. "TRBV11-1*00(1030),TRBV11-3*00(1001)"
# output is the v-segment with the highest score, i.e. "TRBV11-1*00" in this example
get_v = function(v_with_score_single) {
  v_vec = strsplit(v_with_score_single, split = ",")[[1]]
  vs = gsub("\\(.*", "", v_vec)
  scores = gsub(".*\\(([^()]*)\\).*", "\\1", v_vec) %>% as.numeric()
  max_ind = which.max(scores)
  return(as.vector(vs[max_ind]))
}

# input is a vector, a column from a data frame like "AllVHitsWithScore"
# output is a vector of the v-segments with the highest score for each
get_all_vs = function(v_with_score) {
  sapply(v_with_score, function(x) {
    get_v(x)
  }) %>% as.vector()
}
