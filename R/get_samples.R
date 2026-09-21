#' Returns sample names that are available to be imported from a directory
#'
#' @description
#'
#' The \code{get_samples()} takes a directory path and searches the directory for
#' TIRTL-seq files -- files ending in `_TIRTLoutput.tsv(.gz)`, `_pseudobulk_TRA.tsv(.gz)`,
#' or `_pseudobulk_TRB.tsv(.gz)`. You can give a regex pattern to return only
#' samples that match this pattern.
#'
#' @param data_dirs a directory path (or a vector with multiple paths) with TIRTL-seq data files
#' @param pattern a filename pattern
#'
#' @returns
#' A dataframe with a column for the sample name and other columns to indicate
#' if the sample has paired TCR data, TCRalpha pseudobulk data, and/or TCRbeta pseudobulk data.
#'
#' @family data_import
#'
#' @export
#'
#'

get_samples = function(data_dirs, pattern = NULL) {
  df_list = lapply(data_dirs, function(x) get_samples_single_dir(data_dir = x, pattern = pattern))
  return_df = dplyr::bind_rows(df_list)
  ## check for sample name collisions
  dupes = return_df$sample[duplicated(return_df$sample)]
  dupes_char = paste(dupes, collapse = ", ")
  if(length(dupes) > 0) warning(glue::glue("Sample name collision: {dupes_char}"), call. = FALSE)
  return_df = return_df[!duplicated(return_df$sample),]
  return(return_df)
}

#' Internal helper for `get_samples()`
#'
#' @description same as `get_samples()`, but it takes in one directory only
#'
#' @keywords internal
get_samples_single_dir = function(data_dir, pattern = NULL) {
  checkmate::assert_string(data_dir)
  checkmate::assert_directory_exists(data_dir)
  data_dir = sub("/$", "", data_dir) ## remove trailing slash
  all_paired = dir(data_dir, pattern = "_TIRTLoutput\\.tsv") |>
    gsub("_TIRTLoutput\\.tsv.*", "", x = _) |>
    sort()
  all_tra = dir(data_dir, pattern = "_pseudobulk_TRA\\.tsv") |>
    gsub("_pseudobulk_TRA\\.tsv.*", "", x = _) |>
    sort()
  all_trb = dir(data_dir, pattern = "_pseudobulk_TRB\\.tsv") |>
    gsub("_pseudobulk_TRB\\.tsv.*", "", x = _) |>
    sort()
  all_any = c(all_paired, all_tra, all_trb) |> unique() |> sort()
  if(!is.null(pattern)) {
    all_samples = all_any[grepl(pattern, all_any)]
  } else {
    all_samples = all_any
  }
  return_df = tibble::tibble(
    sample = all_samples,
    has_paired = all_samples %in% all_paired,
    has_tra = all_samples %in% all_tra,
    has_trb = all_samples %in% all_trb
    ) |>
    dplyr::mutate(
      has_all = has_paired & has_tra & has_trb,
      char_desc = dplyr::case_when(
        has_paired & has_tra & has_trb ~ "TRA+TRB+paired",
        !has_paired & has_tra & has_trb ~ "TRA+TRB",
        has_paired & has_tra & !has_trb ~ "TRA+paired",
        has_paired & !has_tra & has_trb ~ "TRB+paired",
        has_paired & !has_tra & !has_trb ~ "paired",
        !has_paired & has_tra & !has_trb ~ "TRA",
        !has_paired & !has_tra & has_trb ~ "TRB"
      )
    ) |>
    dplyr::mutate(
      paired_file = sapply(sample, function(x) dir(data_dir, pattern = glue::glue("{x}_TIRTLoutput.tsv.*"), full.names = TRUE)[1]),
      tra_file = sapply(sample, function(x) dir(data_dir, pattern = glue::glue("{x}_pseudobulk_TRA.tsv.*"), full.names = TRUE)[1]),
      trb_file = sapply(sample, function(x) dir(data_dir, pattern = glue::glue("{x}_pseudobulk_TRB.tsv.*"), full.names = TRUE)[1])
    ) ## note: will preferentially choose unzipped ".tsv" files over ".tsv.gz" if both are available due to default sorting in "dir" function

  return(return_df)
}
