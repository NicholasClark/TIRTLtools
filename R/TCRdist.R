#' GPU implementation of TCRdist, a distance/similarity metric for pairs of TCRs
#'
#' @description
#' An efficient, batched version of TCRdist that is compatible with both NVIDIA and Apple Silicon GPUs.
#'
#' @details
#' This function calculates pairwise TCRdist (Dash et al., Nature 2017) for a set of TCRs
#' (or between two sets of TCRs) and returns a sparse output with the TCRdist and indices of all pairs
#' that have TCRdist less than or equal to a desired cutoff (default cutoff is 90).
#'
#' The function uses the `reticulate` package to call a python script that uses `cupy` (NVIDIA GPUs), `mlx` (Apple Silicon GPUs),
#' or `numpy` (no GPU) to calculate TCRdist efficiently.
#'
#' @param tcr1 a data frame with one TCR per row. It must have the columns "va", "vb", "cdr3a", and "cdr3b"
#' @param tcr2 (optional) another data frame of TCRs. If supplied, TCRdist will be calculated
#' for every combination of one TCR from tcr1 and one TCR from tcr2. Otherwise, it will calculate TCRdist
#' for every pair of TCRs in tcr1.
#' @param params (optional) a table of valid parameters for amino acids and va/vb segments.
#' (default is NULL, which uses TIRTLtools::params)
#' @param submat (optional) a substitution matrix with mismatch penalties for each
#' combination of amino acids or va/vb segments (default is NULL, which uses TIRTLtools::submat).
#' @param tcrdist_cutoff (optional) discard all TCRdist values above this cutoff (default is 90).
#' @param chunk_size (optional) The chunk size to use in calculation of TCRdist (default 1000). If set at n,
#' it will calculate pairwise TCRdist for n x n TCRs at once. This may be as high as allowable by GPU memory
#' (in our testing, a chunk_size of 1000 to 5000 provided the fastest runtime and
#' chunk_size of over 7500 resulted in memory errors on some GPUs).
#' @param print_chunk_size (optional) print a line of output for every n TCRs processed (default 1000)
#' @param print_res (optional) print summary of results (default is TRUE)
#' @param only_lower_tri (optional) return one TCRdist value for each pair (like the lower triangle of a symmetric matrix). Default is TRUE.
#' @param return_data (optional) whether to return the output result from the function.
#' With large data it may be desirable to write the result to disk instead. (default is TRUE, returns output)
#' @param write_tsv (optional) write the results to a tab-separated file ".tsv" (default is FALSE, does not write .tsv file)
#'
#' @return
#' A list with entries:
#'
#' \code{$TCRdist_df} - a data frame with three columns: "edge1_0index", "edge2_0index", and "TCRdist".
#' The first two columns contain the indices (0-indexed) of the TCRs for each pair.
#' The last column contains the TCRdist if it is below the cutoff. The output is sparse in that it only contains
#' pairs that have TCRdist <= cutoff.
#'
#' \code{$tcr1} - a data frame of the TCRs supplied to the function. It contains an additional column
#' "tcr_index" with the (0-indexed) index of each TCR.
#'
#' \code{$tcr2} - a similar data frame for tcr2, if it was supplied.
#'
#' @family repertoire_analysis
#' @seealso \code{\link{cluster_tcrs}()}, \code{\link{plot_clusters}()}, and \code{\link{identify_non_functional_seqs}()}
#'
#' @export
#' @examples
#' # example code
#' # data = load_tirtlseq("your_directory/")
#' # df = get_all_tcrs(data, chain="paired", remove_duplicates = TRUE)
#' # out = TCRdist(df, tcrdist_cutoff = 90)

TCRdist = function(
    tcr1,
    tcr2=NULL,
    remove_MAIT = FALSE,
    params = NULL,
    submat = NULL,
    tcrdist_cutoff=90,
    chunk_size=1000,
    print_chunk_size=10,
    print_res = TRUE,
    only_lower_tri = TRUE,
    return_data = TRUE,
    write_to_tsv = FALSE
    ) {
  tcr1 = prep_for_tcrdist(tcr1, params = params, remove_MAIT = remove_MAIT)
  if(!is.null(tcr2)) tcr2 = prep_for_tcrdist(tcr2, params = params, remove_MAIT = remove_MAIT)
  chunk_size = as.integer(chunk_size)
  print_chunk_size = as.integer(print_chunk_size)
  ### load python packages/scripts
  TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = system.file("python", package = "TIRTLtools"), convert = TRUE, delay_load = TRUE)
  pd = reticulate::import("pandas", delay_load = TRUE)
  np = reticulate::import("numpy", delay_load = TRUE)
  ### load substitution matrix and dataframe of parameters
  if(is.null(submat)) submat = TIRTLtools::submat
  if(is.null(params)) params = TIRTLtools::params
  ### convert R objects to python objects
  submat_py = reticulate::r_to_py(submat, convert = TRUE)
  params_py = pd$DataFrame(data = reticulate::r_to_py(params, convert = TRUE))
  tcr1_py = pd$DataFrame(data = reticulate::r_to_py(tcr1, convert = TRUE))
  if(!is.null(tcr2)) {
    tcr2_py = pd$DataFrame(data = reticulate::r_to_py(tcr2, convert = TRUE))
  } else {
    tcr2_py = reticulate::r_to_py(tcr2, convert = TRUE)
  }
  ### call python TCRdist_batch function
  res = TCRdist_gpu$TCRdist_batch(tcr1 = tcr1_py, tcr2 = tcr2_py, submat = submat_py, params_df = params_py, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_chunk_size = print_chunk_size, print_res = print_res, only_lower_tri = only_lower_tri, return_data = return_data, write_to_tsv = write_to_tsv)
  return(res)
}
