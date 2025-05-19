
#' The input 'data' must be a TIRTL-seq data object loaded with the loading function.
#' i.e. it must be a list of length two, with data$meta and data$data, where
#' data$data is a list of samples, each having data for paired, alpha, and beta chains.
add_single_chain_data = function(data) {
  data_tmp = lapply(data$data, function(x) {
    df_out = add_single_chain_data_simple(x)
    return(list(alpha = x$alpha, beta=x$beta, paired = df_out))
  })
  data$data = data_tmp
  data$is_annotated = TRUE
  return(data)
}

#' The input 'data' must be a list with 3 items, "alpha", "beta", and "paired"
#' This looks up readCounts/Fractions of each single chain in the
add_single_chain_data_simple = function(data) {
  dt_pair = data$paired
  dtA = data$alpha
  dtB = data$beta
  colnames(dtA) = paste("alpha_", colnames(dtA), sep = "")
  colnames(dtB) = paste("beta_", colnames(dtB), sep = "")
  colnames(dtA)[1] = "alpha_nuc"
  colnames(dtB)[1] = "beta_nuc"
  dt_join1 = data.table::merge.data.table(x=dt_pair, y=dtA, by = "alpha_nuc", all.x=TRUE, all.y=FALSE)
  dt_join2 = data.table::merge.data.table(x=dt_join1, y=dtB, by = "beta_nuc", all.x=TRUE, all.y=FALSE)
  return(dt_join2)
}
