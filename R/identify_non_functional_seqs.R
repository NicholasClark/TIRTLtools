identify_non_functional_seqs = function(df) {
  if("cdr3a" %in% colnames(df) && "cdr3b" %in% colnames(df)) {
    ### for paired data
    df_out = df %>%
      mutate(
        alpha_has_stop_codon = grepl("\\*", cdr3a), # stop codon
        alpha_has_frameshift = grepl("_", cdr3a), # frame shift
        beta_has_stop_codon = grepl("\\*", cdr3b),
        beta_has_frameshift = grepl("_", cdr3b),
      ) %>%
      mutate(
        alpha_is_functional = !(alpha_has_stop_codon | alpha_has_frameshift),
        beta_is_functional = !(beta_has_stop_codon | beta_has_frameshift)
      ) %>%
      mutate(is_functional = alpha_is_functional & beta_is_functional)
  } else {
    ### for pseudobulk data
    df_out = df %>%
      mutate(
        has_stop_codon = grepl("\\*", aaSeqCDR3), # stop codon
        has_frameshift = grepl("_", aaSeqCDR3), # frame shift
      ) %>%
      mutate(is_functional = !(has_stop_codon & has_frameshift))
  }
  return(df_out)
}
