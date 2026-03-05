#' Load well-level data written to binary format
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function loads well-level data written to binary files (.parquet and .rds)
#' by the \code{\link{write_well_data_to_binary}()} function. The data contains
#' sparse matrices (well x clone) of read counts for TCRalpha and TCRbeta
#' along with metadata data frames for each clone.
#'
#' Data loaded this way can be used as input to the \code{\link{plot_tshell}()} function.
#'
#' @param folder the directory with ".parquet" and ".rds" files written by \code{\link{write_well_data_to_binary}()}
#' @param prefix a prefix with the sample name that will be prepended to the output file names
#' @param lazy whether to "lazy load" the clone metadata files. This can speed up
#' loading because these files are usually very large (default is TRUE).
#'
#' @return
#' This function returns a list with the following elements:
#' - \code{alpha} - a sparse matrix ("dgCMatrix") of TCRalpha read counts in column-oriented format.
#' See \code{\link{Matrix::`dgCMatrix-class`}()}
#' - \code{beta} - a sparse matrix ("dgCMatrix") of TCRbeta read counts in column-oriented format.
#' See \code{\link{Matrix::`dgCMatrix-class`}()}
#' - \code{alpha_meta} - a data frame of clone metadata for TCRalpha clones. Row i relates to column i of \code{alpha}.
#' - \code{beta_meta} - a data frame of clone metadata for TCRbeta clones. Row i relates to column i of \code{beta}.
#' - \code{col_meta} - a data frame of well metadata (i.e. well names)
#' - \code{alpha_meta_file} - the file for the TCRalpha clone metadata
#' - \code{beta_meta_file} - the file for the TCRbeta clone metadata
#'
#' @family well
#'
#' @export
#'


load_well_counts_binary = function(folder, prefix, lazy = TRUE) {
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
                  alpha_meta = alpha_row, beta_meta = beta_row, well_meta = cols,
                  alpha_meta_file = alpha_row_file, beta_meta_file = beta_row_file
  )
  return(out_list)
}
