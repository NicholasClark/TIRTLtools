#' Remove TCRs with nonfunctional CDR3 amino acid sequences
#'
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @param data either a TIRTLseqData object or a data frame with paired TCRs
#' @param verbose whether to print number of TCRs removed
#'
#' @family data_processing
#'
filter_nonfunctional_TCRs = function(data, verbose = TRUE, version = "data.table") {
  is_list = .is.list.only(data)
  if(is_list) {
    for(i in 1:length(data$data)) {
      message(paste("Sample:", names(data$data)[i]))
      data$data[[i]]$paired = .filter_nonfunctional_TCRs_single(df = data$data[[i]]$paired, verbose = verbose, version = version)
    }
  } else {
    data = .filter_nonfunctional_TCRs_single(df = data, verbose = verbose, version = version)
  }
  return(data)
  # df = .identify_non_functional_seqs_single(df)
  # n_tcr_orig = nrow(df)
  # n_nonfunc = sum(!df$is_functional)
  # pct_nonfunc = signif(100*n_nonfunc/n_tcr_orig,2)
  # msg = paste("Removed ", n_nonfunc %>% .add_commas(), " TCRs with non-functional CDR3 amino acid sequences ", "(", pct_nonfunc, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  # if(verbose) message(msg)
  # df = df %>% filter(is_functional) %>%
  #   select(-alpha_has_stop_codon, -alpha_has_frameshift, -beta_has_stop_codon, -beta_has_frameshift,
  #          -alpha_is_functional, -beta_is_functional, -is_functional)
  # return(df)
}


.filter_nonfunctional_TCRs_single = function(df, verbose = TRUE, version = "data.table") {
  cols = colnames(df)
  df = .identify_non_functional_seqs_single(df)
  n_tcr_orig = nrow(df)
  n_nonfunc = sum(!df$is_functional)
  pct_nonfunc = signif(100*n_nonfunc/n_tcr_orig,2)
  msg = paste("Removed ", n_nonfunc %>% .add_commas(), " TCRs with non-functional CDR3 amino acid sequences ", "(", pct_nonfunc, "%) from a total of ", n_tcr_orig %>% .add_commas(), " TCRs.", sep = "")
  if(verbose) message(msg)
  if(version == "data.table") {
    if(!"data.table" %in% class(df)) {
      print("converting to data.table")
      df = as.data.table(df)
    }
    df = df[is_functional==TRUE, ..cols]
  } else {
    print("using dplyr version for testing")
    df = df %>% filter(is_functional) %>%
      select(-alpha_has_stop_codon, -alpha_has_frameshift, -beta_has_stop_codon, -beta_has_frameshift,
             -alpha_is_functional, -beta_is_functional, -is_functional)
  }
  return(df)
}
