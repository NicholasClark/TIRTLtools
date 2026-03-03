load_wells_counts_binary = function(folder, prefix, lazy = TRUE) {
  alpha_mat_file = file.path(folder, paste(prefix, "alpha.rds", sep="_"))
  beta_mat_file = file.path(folder, paste(prefix, "beta.rds", sep="_"))

  alpha_row_file = file.path(folder, paste(prefix, "alpha.parquet", sep="_"))
  beta_row_file = file.path(folder, paste(prefix, "beta.parquet", sep="_"))

  col_file = file.path(folder, paste(prefix, "cols.parquet", sep="_"))

  alpha_mat = readRDS(alpha_mat_file)
  beta_mat = readRDS(beta_mat_file)

  if(lazy) {
    alpha_row = lazy_load_parquet(alpha_row_file)
    beta_row = lazy_load_parquet(beta_row_file)
    cols = arrow::read_parquet(col_file)
  } else {
    alpha_row = arrow::read_parquet(alpha_row_file)
    beta_row = arrow::read_parquet(beta_row_file)
    cols = arrow::read_parquet(col_file)
  }

  out_list = list(alpha = alpha_mat, beta = beta_mat,
                  alpha_meta = alpha_row, beta_meta = beta_row, col_meta = cols,
                  alpha_meta_file = alpha_row_file, beta_meta_file = beta_row_file
  )
  return(out_list)
}
