plot_clusters = function(obj, n_clusters = 10,seed = 1234) {
  largest_clusters = table(obj$df$cluster) %>% sort() %>% tail(n_clusters) %>% as.data.frame() %>% as_tibble() %>% set_colnames(c("cluster", "n")) %>% mutate(cluster = as.character(cluster)) %>% arrange(desc(n))
  ### testing plotting umap or graph of largest clusters
  df_sub = obj$df %>% filter(cluster %in% largest_clusters$cluster)
  idx_keep = unique(df_sub$idx_1index)
  dist_mat = dist_obj_to_matrix(obj, idx_keep)
  dist_mat_noNA = dist_mat
  dist_mat_noNA[is.na(dist_mat_noNA)] = obj$tcrdist_cutoff
  adj_mat = !is.na(dist_mat)
  mode(adj_mat) = "integer"
  type = c("graph", "umap", "heatmap")
  if("umap" %in% type) {
    set.seed(seed)
    um = uwot::umap(as.dist(dist_mat_noNA))
    df_sub$um1 = um[,1]
    df_sub$um2 = um[,2]
    umap_plot = ggplot(df_sub) + geom_point(aes(x=um1, y=um2, color = as.character(cluster)))
  }
  if("graph" %in% type) {
    set.seed(seed)
    gr_sub = igraph::graph_from_adjacency_matrix(adj_sub, mode = "undirected", weighted = NULL)
    graph_plot = make_graph_many_cluster(df_sub, adj_mat, clusters = largest_clusters$cluster)
  }
  if("heatmap" %in% type) {
    set.seed(seed)
    ann_df = df_sub[,"cluster", drop = FALSE] %>% as.data.frame()
    ann_df$cluster = factor(ann_df$cluster, levels = sort(as.integer(largest_clusters$cluster)))
    ann = ComplexHeatmap::columnAnnotation(df = ann_df)
    hm_plot = ComplexHeatmap::Heatmap(adj_mat, top_annotation = ann)
  }
  return(list(umap = umap_plot, graph = graph_plot, heatmap = hm_plot))
}


plot_igraph = function(g, grps, colbar, vertex.size = 4) {
  plot(g, vertex.size = vertex.size, vertex.label = NA)
  legend('topleft', grps, pch=21, col="#777777", pt.bg=colbar, pt.cex=1, cex=.8,
         bty="n", ncol=1)
}

make_graph_many_cluster = function(df, adj_mat, clusters, pdf_name = "", seed = NULL, color_column = "cluster", vertex.size = 4) {
  df_clus = df %>% filter(cluster %in% clusters)
  # dist_df_sub = adj_mat %>% ## get all edges in the cluster
  #   filter(edge1_clone_name %in% df_clus$clone_name,
  #          edge2_clone_name %in% df_clus$clone_name) %>%
  #   select(edge1_clone_name, edge2_clone_name, TCRdist)
  # g = graph_from_data_frame(dist_df_sub, directed = F, vertices = df_clus)
  g = igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = NULL)

  #comps <- components(g)$membership
  #colbar <- rainbow(max(comps)+1)
  #V(g)$color <- colbar[comps+1]

  grp_vec = as.character(df_clus[[color_column]])
  grps = as.character(unique(grp_vec))
  grps = grps[!is.na(grps)]
  n_grps = length(grps)
  colbar = get_colors_12()[1:length(grps)] %>% setNames(grps)
  igraph::V(g)$color = colbar[grp_vec]
  #igraph::V(g)$size <- degree(g)
  if(pdf_name != "") {
    pdf(pdf_name, width = 9, height = 9)
    if(!is.null(seed)) set.seed(seed)
    plot_igraph(g, grps, colbar, vertex.size = vertex.size)
    dev.off()
  }
  if(!is.null(seed)) set.seed(seed)
  plot_igraph(g, grps, colbar, vertex.size = vertex.size)
}
