
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
  n_valid = dim(dist_input)[1]
  sparse_weight_mat_binary = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
                                                  x=dist_df$weight_binary, symmetric = TRUE,
                                                  dims = c(n_valid, n_valid))
  # sparse_tcrdist_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
  #                                           x=dist_df$TCRdist, symmetric = TRUE,
  #                                           dims = c(n_valid, n_valid))
  #
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
  largest_clusters = table(leiden_clust$membership) %>% sort() %>% tail(10) %>% as.data.frame() %>% as_tibble() %>% set_colnames(c("cluster", "n")) %>% arrange(desc(n))
  ### testing plotting umap or graph of top 10 clusters
  df_sub = dist_input %>% filter(cluster %in% largest_clusters$cluster)
  dist_sub = sparse_weight_mat_binary[df_sub$idx_1index, df_sub$idx_1index]
  um = umap::umap(dist_sub)
  df_sub$um1 = um$layout[,1]
  df_sub$um2 = um$layout[,2]
  ggplot(df_sub) + geom_point(aes(x=um1, y=um2, color = as.character(cluster)))

  gr_sub = igraph::graph_from_adjacency_matrix(dist_sub, mode = "undirected", weighted = NULL)

  make_graph_many_cluster(df_sub, dist_sub, clusters = largest_clusters$cluster)

  #ComplexHeatmap::Heatmap(as.matrix(dist_sub))

  # igraph::V(gr_tmp)$size <- igraph::degree(gr_sub)
  # igraph::V(gr_tmp)$source <- df_tmp$source
  # igraph::V(gr_tmp)$name <- df_tmp$clone_name
  # igraph::V(gr_tmp)$epitope_species <- df_tmp$`Epitope species`
  #
  #
  # ggraph(gr_tmp, layout = "stress")  +
  #   geom_edge_link0(aes(edge_linewidth = weight), edge_colour = "grey66") +
  #   geom_node_point(aes(fill = epitope_species, size = size), shape = 21) +
  #   geom_node_text(aes(filter = size > 12, label = name), family = "serif") +
  #   #scale_fill_manual(values = got_palette) +
  #   scale_edge_width(range = c(0.1, 1)) +
  #   scale_size(range = c(1, 6)) +
  #   theme_graph()


}
