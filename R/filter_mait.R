#' Remove MAIT (Mucosal-Associated Invariant T cells) TCRs
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function uses a heuristic to remove TCRs associated with MAIT cells, which
#' may not be desired for some applications. For example, when creating a network of related TCRs
#' with TCRdist, these TCRs often form a very highly dense sub-network that inflates output file size.
#'
#' @details
#' We use a heiristic where we remove TCRs with V-alpha TRAV1-2 and
#' either J-alpha segment TRAJ33, TRAJ12, or TRAJ20 (see Garner et al. 2023 - https://www.nature.com/articles/s41590-023-01575-1).
#'
#' @param df data frame with paired TCRs
#' @param verbose whether to print number of MAIT TCRs removed
#'
#' @returns A data frame with MAIT TCRs removed
#'
#' @family data_processing
#' @export
#'
filter_mait = function(df, verbose = TRUE) {
  # MAIT cells contain TRAV1-2 and TRAJ33/12/20 (https://www.nature.com/articles/s41590-023-01575-1)
  va1 = "TRAV1-2"
  ja1 = c("TRAJ33", "TRAJ12", "TRAJ20")
  va2 = paste(va1, "*01", sep = "") # with alleles
  ja2 = paste(ja1, "*01", sep = "") # with alleles
  n_tcr_orig = nrow(df)
  df = df %>%
    filter( !(va %in% c(va1,va2) & ja %in% c(ja1,ja2)) )
  n_tcr_final = nrow(df)
  n_mait = n_tcr_orig - n_tcr_final
  pct_mait = signif(100*(n_mait/n_tcr_orig), 2)
  msg = paste("Removed ", n_mait %>% .add_commas(), " MAIT TCRs ", "(", pct_mait, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  return(df)
}
