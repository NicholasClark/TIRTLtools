#' Choose a partner manually for an input chain using T-SHELL
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function takes in an object created by \code{\link{plot_tshell}()} and allows
#' a user to choose the partner for the input chain used to create that object.
#'
#' @param res an object created by \code{\link{plot_tshell}()}
#' @param rank the rank of the potential partner chain to choose
#'
#' @return
#' A one-row data frame in the style of the "paired" data frame from pairing output
#' with the input chain and its chosen partner.
#'
#' @family well
#'
#' @export
#'

choose_pair_manual = function(res, rank) {
  df = res$data$df_top_n[rank,]
  chain = res$call$chain
  if(chain == "alpha") {
    df = df %>% dplyr::rename(cdr3b = aaSeqCDR3, beta_nuc = targetSequences, vb = v, jb = j) %>%
      mutate(cdr3a = res$data$input_meta$aaSeqCDR3, alpha_nuc = res$data$input_meta$targetSequences,
             va = res$data$input_meta$v, ja = res$data$input_meta$j)
  } else {
    df = df %>% dplyr::rename(cdr3a = aaSeqCDR3, alpha_nuc = targetSequences, va = v, ja = j) %>%
      mutate(cdr3b = res$data$input_meta$aaSeqCDR3, beta_nuc = res$data$input_meta$targetSequences,
             vb = res$data$input_meta$v, jb = res$data$input_meta$j)
  }
  df = df %>%
    mutate(alpha_beta = paste(alpha_nuc, beta_nuc, sep = "_"), method = "tshell", score = NA) %>%
    dplyr::rename(ts = t, pval = p, pval_adj = p_adj)
  out_cols = c("wi", "wj", "wij", "alpha_nuc", "beta_nuc", "wa", "wb", "alpha_beta",
               "method", "r", "ts", "pval", "pval_adj", "loss_a_frac", "loss_b_frac",
               "score", "cdr3a", "va", "ja", "cdr3b", "vb", "jb")
  df = df[,out_cols]
  return(df)
}
