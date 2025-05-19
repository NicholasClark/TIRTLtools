plot_igraph = function(g, grps, colbar) {
  plot(g, vertex.size = 2, vertex.label = NA)
  legend('topleft', grps, pch=21, col="#777777", pt.bg=colbar, pt.cex=1, cex=.8,
         bty="n", ncol=1)
}

make_graph_many_cluster = function(df, dist_df, clusters, pdf_name = "", seed = NULL, color_column = "cluster") {
  df_clus = df %>% filter(cluster %in% clusters)
  # dist_df_sub = dist_df %>% ## get all edges in the cluster
  #   filter(edge1_clone_name %in% df_clus$clone_name,
  #          edge2_clone_name %in% df_clus$clone_name) %>%
  #   select(edge1_clone_name, edge2_clone_name, TCRdist)
  # g = graph_from_data_frame(dist_df_sub, directed = F, vertices = df_clus)
  g = igraph::graph_from_adjacency_matrix(dist_df, mode = "undirected", weighted = NULL)

  #comps <- components(g)$membership
  #colbar <- rainbow(max(comps)+1)
  #V(g)$color <- colbar[comps+1]

  grp_vec = as.character(df_clus[[color_column]])
  grps = as.character(unique(grp_vec))
  grps = grps[!is.na(grps)]
  n_grps = length(grps)
  colbar = rainbow(n_grps) %>% setNames(grps)
  igraph::V(g)$color = colbar[grp_vec]
  if(pdf_name != "") {
    pdf(pdf_name, width = 9, height = 9)
    if(!is.null(seed)) set.seed(seed)
    plot_igraph(g, grps, colbar)
    dev.off()
  }
  if(!is.null(seed)) set.seed(seed)
  plot_igraph(g, grps, colbar)
}
