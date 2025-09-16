tcrdist_to_matrix = function(tcr_obj) {
  n_vert = dim(tcr_obj$tcr1)[1]
  dist_df = tcr_obj$TCRdist_df %>% mutate %>%
    mutate(edge1_1index = edge1_0index+1L,
           edge2_1index = edge2_0index+1L
    )
  dist_df_to_matrix(dist_df, n_vert, 1:n_vert)
}

dist_df_to_matrix = function(dist_df, n_vertices_all, idx_keep) {
  ## temporarily switch zeros with -1
  dist_df$TCRdist_mod = ifelse(dist_df$TCRdist == 0L, -1L, dist_df$TCRdist)
  ## convert to sparse matrix
  sparse_tcrdist_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
                                            x=dist_df$TCRdist_mod, symmetric = TRUE,
                                            dims = c(n_vertices_all, n_vertices_all))
  tmp = sparse_tcrdist_mat[idx_keep, idx_keep] %>% as.matrix()
  mode(tmp) <- "integer"
  #tmp[tmp==0L] = cutoff
  tmp[tmp==0L] = NA
  tmp[tmp==-1L] = 0
  return(tmp)
}

dist_obj_to_matrix = function(obj, idx_keep) {
  tmp = dist_df_to_matrix(obj$dist_df, n_vertices_all = dim(obj$df)[1], idx_keep = idx_keep)
  return(tmp)
}
