## testing
# load_example_data()
# data1 = SJTRC_minimal$data$cd8_tp1_v2
# data2 = SJTRC_minimal$data$cd8_tp2_v2
# chain = "beta"
# filter_pairs = FALSE
# remove_nonfunctional = FALSE
# log2fc_cutoff = 3
# sem_cutoff = 2.5
# pseudo1 = 1e-6
# pseudo2 = 1e-6
# smooth_sem = c("window", "none")
# window_size = 30
# end_window_size = 5

get_expanded_clones = function( data1,
                                data2,
                                chain = c("beta", "alpha"),
                                filter_pairs = TRUE,
                                remove_nonfunctional = FALSE,
                                log2fc_cutoff = 3,
                                sem_cutoff = 2.5,
                                pseudo1 = 1e-6,
                                pseudo2 = 1e-6,
                                smooth_sem = c("window", "none"),
                                window_size = 30,
                                end_window_size = 5) {
  checkmate::assert_choice(chain, choices = c("alpha", "beta"))
  chain = chain[1]
  smooth_sem = smooth_sem[1]
  other_chain = ifelse(chain == "alpha", "beta", "alpha")

  df_all = .make_sample_vs_sample_df(data1 = data1, data2 = data2, chain = chain, smooth_sem = smooth_sem,
    pseudo1=pseudo1, pseudo2 = pseudo2, sem_cutoff = sem_cutoff, log2fc_cutoff = log2fc_cutoff,
    window_size = window_size, end_window_size = end_window_size
  )
  paired1 = data1$paired_alt %>% mutate(pairs_from = "x")
  paired2 = data2$paired_alt %>% mutate(pairs_from = "y")

  up_df = df_all %>% filter(sign %in% c("up")) %>%
    arrange(desc(abs(log2FC))) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())
  down_df = df_all %>% filter(sign %in% c("down")) %>%
    arrange(desc(abs(log2FC))) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())

  up_out = .merge_with_pairs(single_chain_df = up_df, direction = "up", paired_df1 = paired1, paired_df2 = paired2,
    chain = chain, filter_pairs = filter_pairs, remove_nonfunctional = remove_nonfunctional)
  down_out = .merge_with_pairs(single_chain_df = down_df, direction = "down", paired_df1 = paired1, paired_df2 = paired2,
    chain = chain, filter_pairs = filter_pairs, remove_nonfunctional = remove_nonfunctional)  

  out_list = list(expanded = up_out, contracted = down_out)
  return(out_list)
}

.merge_with_pairs = function(
  single_chain_df, direction, paired_df1, paired_df2,
  chain, filter_pairs, remove_nonfunctional) {
  
  if(chain == "alpha") {
    cdr3_nt_col = "alpha_nuc"
    cdr3_aa_col = "cdr3a"
    v_col = "va"
    j_col = "ja"
  } else {
    cdr3_nt_col = "beta_nuc"
    cdr3_aa_col = "cdr3b"
    v_col = "vb"
    j_col = "jb"
  }
  join_cols = c(cdr3_nt_col, cdr3_aa_col, v_col, j_col)

  df_select = single_chain_df %>% 
    select(
      sign, log2FC, 
      targetSequences, aaSeqCDR3, v, j,
      #avg.x, avg.y, 
      readFraction.x, readFraction.y,
      sem.x, sem.y, n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readCount.x, readCount.y,
      readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y
      ) %>%
    #dplyr::rename(beta_nuc = targetSequences, cdr3b = aaSeqCDR3, vb = v, jb = j)
    dplyr::rename(!!cdr3_nt_col := targetSequences, !!cdr3_aa_col := aaSeqCDR3, !!v_col := v, !!j_col := j)
  
  if(filter_pairs) {
    ## filter paired data frame from each sample down to the best partner for each beta
    paired_df1 = .clean_pairs_single(paired_df1)
    paired_df2 = .clean_pairs_single(paired_df2)
  }
  paired_df1_select = paired_df1 %>% 
    select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from)
  paired_df2_select = paired_df2 %>% 
    select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from)
  
  
  signif_paired_df1 = df_select %>% inner_join(paired_df1_select, by = join_cols)
  signif_paired_df2 = df_select %>% inner_join(paired_df2_select, by = join_cols)
  if(direction == "up") {
    ## if clones are expanding put pairs from second sample at the top (more likely to find pair there)
    signif_paired_all = bind_rows(signif_paired_df2, signif_paired_df1)
  } else {
    ## if clones are contracting put pairs from first sample at the top (more likely to find pair there)
    signif_paired_all = bind_rows(signif_paired_df1, signif_paired_df2)
  }

  signif_paired_all = signif_paired_all %>%
    mutate(comparison_chain = chain) %>%
    filter(!duplicated(alpha_beta)) %>%
    select(
      sign, log2FC,
      va, ja, cdr3a, cdr3b, vb, jb, #avg.x, avg.y,
      alpha_nuc, beta_nuc,
      comparison_chain,
      readFraction.x, readFraction.y,
      sem.x, sem.y, 
      n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readCount.x, readCount.y,
      readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y,
      alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from) %>%
    mutate(across(starts_with("readCount") | starts_with("max_wells") | 
      starts_with("n_wells") | starts_with("rank"), as.integer)) %>%
    arrange(desc(abs(log2FC)))
  
  ## single-chain output
  signif_sc = df_select %>% 
    mutate(across(starts_with("readCount") | starts_with("max_wells") | 
      starts_with("n_wells") | starts_with("rank"), as.integer)) %>%
    mutate(is_paired = !!sym(cdr3_nt_col) %in% signif_paired_all[[cdr3_nt_col]]) %>%
    #filter(!(!!sym(cdr3_nt_col) %in% signif_paired_all[[cdr3_nt_col]])) %>% ## filter out paired single-chains
    select(sign, log2FC, is_paired, everything()) %>%
    arrange(desc(abs(log2FC)))
  if(remove_nonfunctional) {
    signif_paired_all = signif_paired_all %>% filter(is_functional)
    signif_sc = signif_sc %>%
      mutate(has_stop_codon = grepl("\\*", !!sym(cdr3_aa_col)), has_frameshift = grepl("_", !!sym(cdr3_aa_col))) %>%
      mutate(is_functional = !(has_stop_codon | has_frameshift)) %>%
      filter(is_functional)
  }
  out = list(paired = signif_paired_all, single_chain = signif_sc)
  return(out)
}
