#' @title
#' Identify single chains that are paired in the pseudo-bulk data
#'
#' @description
#' For each sample in the dataset, \code{identify_paired()} annotates the alpha
#' and beta pseudobulk data with the number of distinct pairs each chain is a part of
#' in the paired data as well as a TRUE/FALSE column indicating whether the chain
#' is paired with any partners.
#'
#' @param data a TIRTLseq dataset created by \code{\link{load_tirtlseq}()}
#'
#' @return
#' A dataset similar to that created by \code{\link{load_tirtlseq}()}, but
#' with added columns \code{is_paired} and \code{n_paired} in the alpha and beta
#' pseudobulk data frames.
#'
#' \code{is_paired} is TRUE if the chain is found in the
#' paired data. \code{n_paired} is the number of distinct chains that the particular
#' chain is paired with.
#'
#' @seealso \code{\link{load_tirtlseq}()}
#'
#' @export
#' @examples
#' # example code
#' # paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
#' # paired = identify_paired(paired)
#'
identify_paired = function(data) {
  data$data = lapply(data$data, function(x) identify_paired_single(x))
  return(data)
}


identify_paired_single = function(data) {
  paired_tmp = data$paired %>% remove_dupes_paired()
  counts_alpha = paired_tmp[, .N, by = alpha_nuc] %>% dplyr::rename(n_paired = N, targetSequences = alpha_nuc)
  counts_beta = paired_tmp[, .N, by = beta_nuc] %>% dplyr::rename(n_paired = N, targetSequences = beta_nuc)
  data$alpha = left_join(data$alpha, counts_alpha) %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0)
  data$beta = left_join(data$beta, counts_beta) %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0)
  return(data)
}
