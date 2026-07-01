
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
  merge_df1 = data1$paired_alt
  merge_df2 = data2$paired_alt

  up_df = df_all %>% filter(sign %in% c("up")) %>%
    arrange(desc(log2FC)) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())
  down_df = df_all %>% filter(sign %in% c("down")) %>%
    arrange(log2FC) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())
  if(chain == "beta") {
    up_df_select = up_df %>% select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, 
      #avg.x, avg.y, 
      sem.x, sem.y, n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readFraction.x, readFraction.y, readCount.x, readCount.y,
      readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y
      ) %>%
      dplyr::rename(beta_nuc = targetSequences, cdr3b = aaSeqCDR3, vb = v, jb = j)
    down_df_select = down_df %>% select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, 
      #avg.x, avg.y, 
      sem.x, sem.y, n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readFraction.x, readFraction.y, readCount.x, readCount.y) %>%
      dplyr::rename(beta_nuc = targetSequences, cdr3b = aaSeqCDR3, vb = v, jb = j)

    paired1 = data1$paired_alt
    paired2 = data2$paired_alt
    if(filter_pairs) {
      ## filter paired data frame from each sample down to the best partner for each beta
      paired1 = .clean_pairs_single(paired1)
      paired2 = .clean_pairs_single(paired2)
    }
    paired1_select = paired1 %>% select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional)
    paired2_select = paired2 %>% select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional)
    
    
    up_paired_df1 = up_df_select %>% inner_join(paired1_select, by = c("beta_nuc", "cdr3b", "vb", "jb"))
    up_paired_df2 = up_df_select %>% inner_join(paired2_select, by = c("beta_nuc", "cdr3b", "vb", "jb"))
    up_paired_all = bind_rows(up_paired_df2, up_paired_df1) %>%
      mutate(comparison_chain = chain) %>%
      #filter(!duplicated(paste(beta_nuc, cdr3b, vb, jb))) %>%
      select(va, ja, cdr3a, cdr3b, vb, jb, comparison_chain, sign, log2FC, #avg.x, avg.y,
        readFraction.x, readFraction.y, sem.x, sem.y, 
        n_wells.x, n_wells.y, max_wells.x, max_wells.y,
        rank.x, rank.y, readCount.x, readCount.y,
        readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y,
        alpha_nuc, beta_nuc, alpha_beta, 
        is_functional, alpha_is_functional, beta_is_functional) %>%
      mutate(n_wells.x = as.integer(n_wells.x), n_wells.y = as.integer(n_wells.y),
        max_wells.x = as.integer(max_wells.x), max_wells.y = as.integer(max_wells.y),
        rank.x = as.integer(rank.x), rank.y = as.integer(rank.y)
        ) %>%
      arrange(desc(log2FC))
    
    up_unpaired = up_df_select %>% filter(!beta_nuc %in% up_paired_all$beta_nuc) %>% 
      mutate(n_wells.x = as.integer(n_wells.x), n_wells.y = as.integer(n_wells.y),
        max_wells.x = as.integer(max_wells.x), max_wells.y = as.integer(max_wells.y),
        rank.x = as.integer(rank.x), rank.y = as.integer(rank.y)
        ) %>%
      arrange(desc(log2FC))
    if(remove_nonfunctional) {
      up_paired_all = up_paired_all %>% filter(is_functional)
      up_unpaired = up_unpaired %>%
        mutate(has_stop_codon = grepl("\\*", cdr3b), has_frameshift = grepl("_", cdr3b)) %>%
        mutate(is_functional = !(has_stop_codon | has_frameshift)) %>%
        filter(is_functional)
    }
    down_paired_df1 = data1$paired_alt %>% inner_join(down_df_select, by = c("beta_nuc", "cdr3b", "vb", "jb"))
    down_paired_df2 = data2$paired_alt %>% inner_join(down_df_select, by = c("beta_nuc", "cdr3b", "vb", "jb"))
    
  } else {
    up_df_select = up_df %>% select(targetSequences, aaSeqCDR3, v, j, sign, log2FC) %>%
      dplyr::rename(alpha_nuc = targetSequences, cdr3a = aaSeqCDR3, va = v, ja = j) %>%
      mutate(match_col = paste(alpha_nuc, cdr3a, va, ja))
  }
  
  
  

  out_list = list(expanded = list(single_chain = up_df, paired_data1 = up_paired_df1, paired_data2 = up_paired_df2),
                  contracted = list(single_chain = down_df, paired_data1 = down_paired_df1, paired_data2 = down_paired_df2)
                )
  return(out_list)
}