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
identify_paired = function(data, verbose = TRUE) {
  data$data = lapply(1:length(data$data), function(i) {
    if(verbose) {
      msg = paste("Annotating data with pairing status by MAD-HYPE and T-SHELL algorithms for sample", i) %>% .prepend_newline()
      cat(msg)
    }
    x=data$data[[i]]
    .annotate_paired_single(x)
    }) %>% set_names(names(data$data))
  return(data)
}

## input is data for one sample/experiment
.get_pair_stats_single = function(data) {
  if(!"is_paired" %in% colnames(data$beta)) {
    data = .annotate_paired_single(data)
  }
  paired_df = data$paired_alt$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq")) %>%
    mutate(chain = "paired")
  alpha_df = data$alpha$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq")) %>%
    mutate(chain = "alpha")
  beta_df = data$beta$paired_status %>% table() %>%
    as.data.frame.table() %>% as_tibble() %>%
    set_colnames(c("category", "Freq"))  %>%
    mutate(chain = "beta")
  all_df = bind_rows(list(paired_df, alpha_df, beta_df))
  return(all_df)
}


.annotate_paired_single = function(data) {
  ll = .identify_paired_single(data)
  data$alpha = ll$alpha
  data$beta = ll$beta
  data$paired_alt = ll$paired_alt ## paired data frame w/ no duplicates
  return(data)
}

# for input of a single sample/experiment, returns alpha, beta, and
# paired (with no duplicates) data frames with methods (madhype and tshell) annotated
.identify_paired_single = function(data) {
  paired_tmp = data$paired %>% remove_dupes_paired()
  paired_tmp_madhype = data$paired %>% filter(method == "madhype")
  paired_tmp_tshell = data$paired %>% filter(method == "tshell")

  counts_alpha_madhype = paired_tmp_madhype[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired_madhype = N, targetSequences = alpha_nuc)
  counts_beta_madhype = paired_tmp_madhype[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired_madhype = N, targetSequences = beta_nuc)


  counts_alpha_tshell = paired_tmp_tshell[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired_tshell = N, targetSequences = alpha_nuc)
  counts_beta_tshell = paired_tmp_tshell[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired_tshell = N, targetSequences = beta_nuc)

  counts_alpha = paired_tmp[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired = N, targetSequences = alpha_nuc)
  counts_beta = paired_tmp[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired = N, targetSequences = beta_nuc)

  tmp_alpha = data$alpha %>%
    left_join(counts_alpha_madhype, by = "targetSequences") %>%
    left_join(counts_alpha_tshell, by = "targetSequences") %>%
    left_join(counts_alpha, by = "targetSequences") %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0) %>%
    mutate(n_paired_tshell = ifelse(is.na(n_paired_tshell), 0, n_paired_tshell)) %>%
    mutate(is_paired_tshell = n_paired_tshell != 0) %>%
    mutate(n_paired_madhype = ifelse(is.na(n_paired_madhype), 0, n_paired_madhype)) %>%
    mutate(is_paired_madhype = n_paired_madhype != 0) %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    ))

  tmp_beta = data$beta %>%
    left_join(counts_beta_madhype, by = "targetSequences") %>%
    left_join(counts_beta_tshell, by = "targetSequences") %>%
    left_join(counts_beta, by = "targetSequences") %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0) %>%
    mutate(n_paired_tshell = ifelse(is.na(n_paired_tshell), 0, n_paired_tshell)) %>%
    mutate(is_paired_tshell = n_paired_tshell != 0) %>%
    mutate(n_paired_madhype = ifelse(is.na(n_paired_madhype), 0, n_paired_madhype)) %>%
    mutate(is_paired_madhype = n_paired_madhype != 0)  %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    ))

  paired_tmp = paired_tmp %>%
    mutate(method = NULL) %>%
    mutate(is_paired_madhype = alpha_beta %in% paired_tmp_madhype$alpha_beta) %>%
    mutate(is_paired_tshell = alpha_beta %in% paired_tmp_tshell$alpha_beta) %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    )) ## paired data frame with no duplicates, labeled by method

  return(list(alpha = tmp_alpha, beta = tmp_beta, paired_alt = paired_tmp))
}
