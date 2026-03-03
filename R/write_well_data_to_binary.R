write_well_data_to_binary = function(
    folder_in,
    folder_out,
    prefix,
    wells = get_well_subset(1:16,1:24),
    in_file_type = c(".tsv", ".parquet"),
    out_file_type = c(".rds", ".h5"),
    well_pos=3,
    parallel = FALSE,
    nproc = data.table::getDTthreads(),
    columns = NULL,
    max_files = Inf,
    periods_to_underscores = TRUE,
    to_sparse_matrix = TRUE
) {

  .create_folder(folder_out)
  df_alpha = load_well_data_to_df_single(
    folder_in,
    wells = wells,
    chain = "alpha",
    file_type = in_file_type,
    well_pos = well_pos,
    parallel = parallel,
    nproc = nproc,
    columns = columns,
    max_files = max_files,
    periods_to_underscores = periods_to_underscores
  )
  df_beta = load_well_data_to_df_single(
    folder_in,
    wells = wells,
    chain = "beta",
    file_type = in_file_type,
    well_pos = well_pos,
    parallel = parallel,
    nproc = nproc,
    columns = columns,
    max_files = max_files,
    periods_to_underscores = periods_to_underscores
  )
  if(to_sparse_matrix) {
    ll_alpha = well_data_df_to_sparse(df_alpha, cols = wells, transpose = T)
    ll_beta = well_data_df_to_sparse(df_beta, cols = wells, transpose = T)
    check = all.equal(ll_alpha$col_data, ll_beta$col_data)
    if(!check) stop("well ordering for alpha and beta doesn't match")
    if(!out_file_type %in% c(".rds", ".h5")) stop("'out_file_type' must be either '.rds' or '.h5'")
    if(out_file_type == ".h5") {
      HDF5Array::writeTENxMatrix(
        ll_alpha$matrix,
        filepath = file.path(folder_out, paste(prefix,"alpha.h5",sep="_")),
        group = "matrix")
      HDF5Array::writeTENxMatrix(
        ll_beta$matrix,
        filepath = file.path(folder_out, paste(prefix,"beta.h5",sep="_")),
        group = "matrix")
    } else if(out_file_type == ".rds") {
      saveRDS(ll_alpha$matrix, file = file.path(folder_out, paste(prefix,"alpha.rds",sep="_")))
      saveRDS(ll_beta$matrix, file = file.path(folder_out, paste(prefix,"beta.rds",sep="_")))
    }
    arrow::write_parquet(
      ll_alpha$row_data,
      file.path(folder_out, paste(prefix,"alpha.parquet",sep="_")))
    arrow::write_parquet(
      ll_beta$row_data,
      file.path(folder_out, paste(prefix,"beta.parquet",sep="_")))
    arrow::write_parquet(
      ll_beta$col_data,
      file.path(folder_out, paste(prefix,"cols.parquet",sep="_")))
  } else {
    arrow::write_parquet(
      df_alpha,
      file.path(folder_out, paste(prefix,"alpha.parquet",sep="_")))
    arrow::write_parquet(
      df_beta,
      file.path(folder_out, paste(prefix,"beta.parquet",sep="_")))
  }
  return(invisible(NULL))
}
