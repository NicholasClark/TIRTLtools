

### Annotate a clustering by the most common antigen species, genes, etc.

## inputs -- df - a data frame with a column noting the cluster (default column name: "cluster"), a "source" column with "observed" for non-annotated rows
annotate_clusters = function(df, cluster_col = "cluster", min_cluster_size = 5, min_annotations = 2, min_obs = 1, db_cols = c("antigen.species", "antigen.gene", "antigen.epitope", "mhc.a", "mhc.b", "mhc.class")) {
  df_summ = df %>% group_by(!!sym(cluster_col)) %>%
    summarize(cluster_size = n(), n_obs = sum(source == "observed"), n_annotated = sum(source != "observed"))



  # clusters_singletons = df_summ %>% filter(cluster_size == 1) %>% extract2(cluster_col)
  # df_singletons = df %>% filter(!!sym(cluster_col) %in% cluster_singletons)
  # clusters_unannotated = df_summ %>% filter(n_annotated == 0) %>% extract2(cluster_col)
  # df_unannotated = df %>% filter(!!sym(cluster_col) %in% clusters_unannotated)
  # clusters_no_obs = df_summ %>% filter(n_obs == 0) %>% extract2(cluster_col)
  # df_no_obs = df %>% filter(!!sym(cluster_col) %in% clusters_no_obs)
  # df_keep = df[!df[[cluster_col]] %in% unique(c(clusters_singletons, clusters_unannotated, clusters_no_obs)),]
  # df_send = df_keep %>% filter(cluster_size >= min_cluster_size, n_annotated >= min_annotations, n_obs >= min_obs)

  df_summ_keep = df_summ %>% filter(cluster_size >= min_cluster_size, n_annotated >= min_annotations, n_obs >= min_obs)
  df_send = df %>% filter(cluster %in% df_summ_keep$cluster)
  res_df = get_all_annotations(df_send, clusters = unique(df_send[[cluster_col]]), cluster_col = cluster_col, db_cols = db_cols)
  return(res_df)
}

get_all_annotations = function(df, clusters, cluster_col = "cluster", db_cols = c("antigen.species", "antigen.gene", "antigen.epitope", "mhc.a", "mhc.b", "mhc.class")) {
  summ_df = lapply(clusters, function(clus) {
    df_clus = df %>% filter(!!sym(cluster_col) == clus)
    get_annotations(df_clus, cluster_col, db_cols) %>% mutate(!!sym(cluster_col):= clus)
  }) %>% bind_rows() %>%
    select(!!sym(cluster_col), everything())
  return(summ_df)
}

### Find most common annotations for many columns
get_annotations = function(df_clus, cluster_col = "cluster",
                           db_cols = c("antigen.species", "antigen.gene", "antigen.epitope", "mhc.a", "mhc.b", "mhc.class")) {
  ll = lapply(1:length(db_cols), function(i) summarize_db(df_clus, col = db_cols[i], cluster_col = cluster_col, include_columns = i==1))
  df_out = bind_cols(ll)
  return(df_out)
}


### Find the most common annotations in each cluster
summarize_db = function(df_clus, col, cluster_col = "cluster", include_columns = TRUE) {
  n_total = nrow(df_clus)
  n_obs = sum(df_clus$source == "observed")
  n_annotated = sum(df_clus$source != "observed")
  df_clus_ann = df_clus %>% filter(source != "observed") %>%
    mutate(full_tcr = paste(va, ja, cdr3a, cdr3b, vb, jb))
  tbl = table(df_clus_ann[[col]]) %>% sort() %>% rev()
  df_obs = df_clus %>% filter(source == "observed")
  n_obs_unique = length(unique(paste(df_obs$alpha_nuc, df_obs$beta_nuc, df_obs$va, df_obs$vb, df_obs$ja, df_obs$jb)))
  n_annotated_unique = length(unique(paste(df_clus_ann$alpha_nuc, df_clus_ann$beta_nuc, df_clus_ann$va, df_clus_ann$vb, df_clus_ann$ja, df_clus_ann$jb)))

  tmp = tibble(cluster = unique(df_clus[[cluster_col]]), n_total = n_total, n_obs = n_obs,
               n_obs_unique = n_obs_unique, n_annotated = n_annotated, n_annotated_unique = n_annotated_unique)

  n_start = ncol(tmp)+1

  tmp[[paste("most_common", col, sep = "_")]] = names(tbl)[1]
  tmp[[paste("num_most_common", col, sep = "_")]] = tbl[1]
  tmp[[paste("percent_most_common", col, sep = "_")]] = 100*(tbl[1]/sum(tbl))
  tmp[[paste("num_other", col, sep = "_")]] = sum(tbl[-1])
  if(!include_columns) tmp = tmp[,n_start:ncol(tmp)]
  return(tmp)
}



