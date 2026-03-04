#' Convert individual well .tsv files to faster loading binary format
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function takes in a folder of .tsv files with TCRalpha and TCRbeta
#' read counts of individual wells and converts the data from a subset of these wells
#' (or all wells) to sparse matrices (well x clone)
#' of read counts along with metadata data frames for each clone.
#'
#' This output can then be quickly loaded with \code{\link{load_well_counts_binary}()}
#' and used in the future as input to
#' the \code{\link{plot_tshell}()} function.
#'
#' @param folder_in the directory with ".tsv" files with read counts for each well
#' @param folder_out the directory to write the output data to (if this does not exist,
#' it will be created)
#' @param prefix a prefix with the sample name that will be prepended to the output file names
#' @param wells a vector of the wells corresponding to the sample (default is all wells on the
#' 384-well plate)
#' @param in_file_type the type of files in \code{folder_in} (default is ".tsv")
#' @param out_file_type the type of output file for the sparse read count matrices
#' for TCRalpha and TCRbeta (default is ".rds")
#' @param well_pos the position of the well ID (e.g. "B5") in the file names
#' (when separating by underscores). For example, files
#' named "<lab>_<project>_<well_id>_TCRalpha.tsv" would use well_pos=3. (default is 3)
#' @param parallel whether to use multiple processors when loading the data (default is FALSE)
#' @param nproc number of processors to use when parallel is TRUE
#' @param columns columns of files to read (default is all columns)
#' @param max_files (for testing purposes) the maximum number of files to load (default is all files)
#' @param periods_to_underscores whether to convert periods to underscores in
#' well files names for determining the well name from \code{well_pos}. Default is TRUE to
#' allow for compatibility with some pre-existing datasets.
#' @param to_sparse_matrix whether to write the read counts to a sparse matrix
#' vs. a long data frame (default is TRUE)
#'
#' @return
#' The function writes the data to the output folder and then returns NULL.
#'
#' @family well
#'
#' @export
#'

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
