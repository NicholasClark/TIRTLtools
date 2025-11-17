#' Remove excess pairs for individual single chains
#'
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function filters the paired TCR data, keeping at most one beta-chain partner for each
#' individual alpha-chain and at most two alpha-chain partners for each individual beta-chain.
#'
#' @details Excess partners for an individual chain are often sequencing errors of the true partner.
#' The sequencing error-derived chains often share the same V/J segment as the true partner, but are
#' found at much lower read fractions.
#'
#' As a heuristic, to mitigate this phenomenon,
#' for each unique beta chain we group the partner alpha chains by their V/J segments and keep only
#' the alpha chains with the highest read fraction for each group. Out of the remaining alpha chain
#' partners, we keep up to two chains (can be changed with n_max_alpha) with the highest read fractions.
#'
#' We go through a similar process with each unique alpha chain, grouping partner beta
#' chains by their V/J segments. However, we keep only one beta chain (can be changed with n_max_beta).
#'
#' @param data a TIRTLseqData object
#' @param n_max_alpha (optional) the maximum number of alpha chains allowed
#' paired with an individual beta chain (default 2)
#' @param n_max_beta (optional) the maximum number of beta chains allowed
#' paired with an individual alpha chain (default 1)
#' @param verbose (optional) whether to print progress of the function (default is TRUE).
#' @param prefer_functional (optional) if TRUE, prefer functional to non-functional chains even if they are less frequent
#' @param version (optional) For testing purposes -- if "fast" (default) use fast data.table version of function, otherwise use dplyr version.
#'
#' @family data_processing
#'
clean_pairs = function(data, n_max_alpha = 2, n_max_beta = 1, verbose = TRUE, prefer_functional = TRUE, version = "fast") {
  if(is.null(data$is_annotated)) {
    data = add_single_chain_data(data, verbose = verbose)
  }
  data_tmp = lapply(1:length(data$data), function(i) {
    x = data$data[[i]]
    if(verbose) {
      msg = paste("Cleaning paired dataframe for sample", i) %>% .add_newline()
      cat(msg)
    }
    paired_out = .clean_pairs_single(x$paired, n_max_alpha = n_max_alpha, n_max_beta = n_max_beta,
                                     prefer_functional = prefer_functional, version = version)
    if(!is.null(x$paired_alt)) {
      paired_alt_out = .clean_pairs_single(x$paired_alt, n_max_alpha = n_max_alpha, n_max_beta = n_max_beta,
                                           prefer_functional = prefer_functional, version = version)
      return(list(alpha = x$alpha, beta=x$beta, paired = paired_out, paired_alt = paired_alt_out))
    } else {
      return(list(alpha = x$alpha, beta=x$beta, paired = paired_out))
    }
  }) %>% set_names(names(data$data))
  data$data = data_tmp
  data$pairs_cleaned = TRUE
  return(data)
}



# input is one data frame of paired TCRs
# output is a "cleaned" data frame where each alpha chain has its matches pruned to
# at most one beta chain and each beta chain has its matches pruned to at most two alpha chains.
.clean_pairs_single = function(df, n_max_alpha = 2, n_max_beta = 1, prefer_functional = TRUE, version = "fast") {
  #df_orig = df
  cols = colnames(df)
  df = df %>% remove_duplicates() %>% .identify_non_functional_seqs_single()
  if(!"data.table" %in% class(df)) {
    print("converting to data.table")
    df = as.data.table(df)
  }
  if(prefer_functional) {
    ## Fast version with data.table
    if(version == "fast") {
      cleaned_df = df[
        order(beta_nuc, vb, jb, va, ja, -is_functional, -alpha_readFraction),
        head(.SD, 1),
        by = .(beta_nuc, vb, jb, va, ja)
      ][
        order(beta_nuc, vb, jb, -is_functional, -alpha_readFraction),
        head(.SD, n_max_alpha),
        by = .(beta_nuc, vb, jb)
      ][
        order(alpha_nuc, va, ja, vb, jb, -is_functional, -beta_readFraction),
        head(.SD, 1),
        by = .(alpha_nuc, va, ja, vb, jb)
      ][
        order(alpha_nuc, va, ja, -is_functional, -beta_readFraction),
        head(.SD, n_max_beta),
        by = .(alpha_nuc, va, ja)
      ]
    } else {
      # Slow version with dplyr for checking
      print("using dplyr version (for testing)")
      cleaned_df = df %>%
        group_by(beta_nuc, vb, jb, va, ja) %>%
        ## if multiple alphas with same va + ja are paired with same beta, assume most frequent is real and rest are sequencing errors
        arrange(desc(is_functional), desc(alpha_readFraction), .by_group = TRUE) %>%
        slice_head(n=1) %>%
        ungroup() %>%
        group_by(beta_nuc, vb, jb) %>%
        arrange(desc(is_functional), desc(alpha_readFraction), .by_group = TRUE) %>%
        slice_head(n=n_max_alpha) %>%
        ungroup() %>%
        ## if multiple betas with same vb + jb are paired with same alpha, assume most frequent is real and rest are sequencing errors
        group_by(alpha_nuc, va, ja, vb, jb) %>%
        arrange(desc(is_functional), desc(beta_readFraction), .by_group = TRUE) %>%
        slice_head(n=1) %>%
        ungroup() %>%
        group_by(alpha_nuc, va, ja) %>%
        arrange(desc(is_functional), desc(beta_readFraction), .by_group = TRUE) %>%
        slice_head(n=n_max_beta) %>%
        ungroup()
    }
  } else {
    if(version == "fast") {
      ## Fast version with data.table
      cleaned_df = df[
        order(beta_nuc, vb, jb, va, ja, -alpha_readFraction),
        head(.SD, 1),
        by = .(beta_nuc, vb, jb, va, ja)
      ][
        order(beta_nuc, vb, jb, -alpha_readFraction),
        head(.SD, n_max_alpha),
        by = .(beta_nuc, vb, jb)
      ][
        order(alpha_nuc, va, ja, vb, jb, -beta_readFraction),
        head(.SD, 1),
        by = .(alpha_nuc, va, ja, vb, jb)
      ][
        order(alpha_nuc, va, ja, -beta_readFraction),
        head(.SD, n_max_beta),
        by = .(alpha_nuc, va, ja)
      ]
    } else {
      print("using dplyr version (for testing)")
      ## Slow version with dplyr for checking
      cleaned_df = df %>%
        group_by(beta_nuc, vb, jb, va, ja) %>%
        ## if multiple alphas with same va + ja are paired with same beta, assume most frequent is real and rest are sequencing errors
        slice_max(alpha_readFraction, n=1) %>% ## take top alpha for each unique beta + va + ja
        ungroup() %>%
        group_by(beta_nuc, vb, jb) %>%
        slice_max(alpha_readFraction, n=n_max_alpha) %>% ## take at most 2 of those alphas for each beta
        ungroup() %>%
        ## if multiple betas with same vb + jb are paired with same alpha, assume most frequent is real and rest are sequencing errors
        group_by(alpha_nuc, va, ja, vb, jb) %>%
        slice_max(beta_readFraction, n=1) %>% ## take top beta for each unique alpha + vb + jb
        ungroup() %>%
        group_by(alpha_nuc, va, ja) %>%
        slice_max(beta_readFraction, n=n_max_beta) %>% ## take top beta for each alpha
        ungroup()
    }
  }
  if(!"data.table" %in% class(cleaned_df)) {
    print("converting to data.table")
    cleaned_df = as.data.table(cleaned_df)
  }
  return(cleaned_df[,..cols])
}

## input is one data frame of paired TCRs
.count_pairs_single = function(df) {
  ### one alpha with many betas
  multi_alpha = df %>%
    group_by(alpha_nuc, cdr3a, va, ja) %>%
    summarize(n_pairs_alpha = n()) %>%
    ungroup() %>%
    arrange(desc(n_pairs_alpha))

  df_alpha = df %>%
    left_join(multi_alpha) %>%
    arrange(desc(n_pairs_alpha), desc(cdr3a), desc(beta_readCount)) %>%
    select(cdr3b, vb, jb, cdr3a, va, ja, alpha_readFraction, beta_readFraction, alpha_readCount, beta_readCount, wa, wb, n_pairs_alpha, beta_nuc, alpha_nuc, everything())

  ### one beta with many alphas
  multi_beta = df %>%
    group_by(beta_nuc, cdr3b, vb, jb) %>%
    summarize(n_pairs_beta = n()) %>%
    ungroup() %>%
    arrange(desc(n_pairs_beta))

  df_beta = df %>%
    left_join(multi_beta) %>%
    arrange(desc(n_pairs_beta), desc(cdr3b), desc(alpha_readCount)) %>%
    select(cdr3b, vb, jb, cdr3a, va, ja, alpha_readFraction, beta_readFraction, alpha_readCount, beta_readCount, wa, wb, n_pairs_beta, beta_nuc, alpha_nuc, everything())

  df_ret = df %>%
    left_join(multi_beta) %>%
    left_join(multi_alpha) %>%
    arrange(desc(n_pairs_beta), desc(cdr3b), desc(alpha_readCount)) %>%
    select(cdr3b, vb, jb, cdr3a, va, ja, alpha_readFraction, beta_readFraction, alpha_readCount, beta_readCount, wa, wb, n_pairs_alpha, n_pairs_beta, beta_nuc, alpha_nuc, everything())
}
