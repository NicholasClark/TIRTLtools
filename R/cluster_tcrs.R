
cluster_tcrs = function(data, cluster.type = "leiden", tcrdist_cutoff = 90, resolution = 0.1) {
  chain = "paired"
  df_all = get_all_tcrs(data, chain, remove_duplicates = TRUE) ## get all tcrs in one data frame
  #df_all = df_all %>% filter(is_functional)
  dist = TCRdist(df_all, tcrdist_cutoff = tcrdist_cutoff)
  dist_df = dist$TCRdist_df
  dist_input = dist$tcr1

  dist_df = dist_df %>% mutate(
    edge1_1index = edge1_0index + 1,
    edge2_1index = edge2_0index + 1
  )
  dist_df$weight_binary = 1
  #dist_df$TCRdist_mod = ifelse(dist_df$TCRdist_mod == 0, -1, dist_df$TCRdist_mod)
  n_valid = dim(dist_input)[1]
  sparse_weight_mat_binary = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
                                                  x=dist_df$weight_binary, symmetric = TRUE,
                                                  dims = c(n_valid, n_valid))
  # sparse_tcrdist_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
  #                                           x=dist_df$TCRdist_mod, symmetric = TRUE,
  #                                           dims = c(n_valid, n_valid))

  # sparse_weight_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
  #                                          x=dist_df$weight, symmetric = TRUE,
  #                                          dims = c(n_valid, n_valid) )

  gr_binary = igraph::graph_from_adjacency_matrix(sparse_weight_mat_binary, mode = "undirected", weighted = NULL)
  leiden_clust = igraph::cluster_leiden(gr_binary, resolution = resolution)
  dist_input$idx_1index = 1:dim(dist_input)[1]
  dist_input$cluster = leiden_clust$membership ### assign clusters to TCR sequences
  tab = table(leiden_clust$membership)
  single = sum(tab == 1)
  g2 = sum(tab >= 2)
  g10 = sum(tab >= 10)
  g50 = sum(tab >= 50)
  g100 = sum(tab >= 100)
  msg1 = paste("Out of ", dim(dist_input)[1], " valid TCRs, ", g2, " clusters detected and ", single, " singleton TCRs.", sep = "") %>% add_newline()
  msg2 = paste(g10, " clusters of size >= 10, ", g50, " clusters of size >= 50, ", g100, " clusters of size >=100.", sep = "") %>% add_newline()
  cat(msg1); cat(msg2)
  out = list(
    df = dist_input,
    dist_df = dist_df,
    #sparse_tcrdist_mat = sparse_tcrdist_mat, ### returning this as a sparse matrix is problematic because missing entries will be seen as TCRdist = 0 instead of TCRdist > cutoff.
    sparse_adj_mat = sparse_weight_mat_binary,
    graph_adj = gr_binary,
    tcrdist_cutoff = tcrdist_cutoff,
    resolution = resolution
    )
  return(out)
}
