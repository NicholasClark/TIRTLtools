### plot network of multi-pairing

plot_multipairing = function(df) {
  df = remove_duplicates(df) %>%
    mutate(alpha_nuc = paste("alpha_",alpha_nuc, sep = ""),
           beta_nuc = paste("beta_", beta_nuc, sep = ""))
  edges = df[,c("alpha_nuc", "beta_nuc")]
  mp_alpha = table(df$alpha_nuc) %>% sort() %>% rev() %>% as.data.frame()
  mp_beta = table(df$beta_nuc) %>% sort() %>% rev() %>% as.data.frame()
  mp_alpha_sub = mp_alpha %>% filter(Freq > 5)
  mp_beta_sub = mp_beta %>% filter(Freq > 5)
  edges_sub = edges %>% filter(alpha_nuc %in% mp_alpha_sub$Var1 | beta_nuc %in% mp_beta_sub$Var1)
  g = igraph::graph_from_data_frame(edges_sub, directed = FALSE)
  V(g)$type = grepl("^alpha", V(g)$name)
  V(g)$chain = ifelse(V(g)$type, "alpha", "beta")
  plot(g, vertex.label=NA, vertex.size=3, vertex.color = ifelse(V(g)$type, "red", "blue") )
  #comp <- components(g)
  #keep_comps <- which(comp$csize >= 3)
  #keep_vertices <- V(g)[comp$membership %in% keep_comps]
  #g_sub <- induced_subgraph(g, vids = keep_vertices)

  # edge_list <- as_data_frame(g, what = "edges")
  # node_list <- as_data_frame(g, what = "vertices")
  # node_list = node_list %>% mutate(chain = ifelse(grepl("^alpha", name), "alpha", "beta"))
  # write_tsv(edge_list, file = "~/edges_pairing_test.tsv")
  # write_tsv(node_list, file = "~/nodes_pairing_test.tsv")

  i=2
  df_ia = df %>% filter(alpha_nuc %in% mp_alpha_sub$Var1[i])
  df_ib = df %>% filter(beta_nuc %in% mp_beta_sub$Var1[i])
  aph63 = load_wells_counts_h5("mvp218_aph63_cd8_ud21", folder = "~/git/multipairing_exploration/data")
  #comp <- igraph::components(g)


}
