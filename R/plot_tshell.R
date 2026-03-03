plot_tshell = function(
    well_data,
    nuc_seq,
    chain,
    n_plot = 9,
    interactive = FALSE,
    normalize_counts = TRUE,
    wij_cutoff = 0.5
) {
  if(!chain %in% c("alpha", "beta")) stop("'chain' needs to be 'alpha' or 'beta'")

  alpha_mat = well_data$alpha
  beta_mat = well_data$beta

  alpha_mat_rf = normalize_rows(alpha_mat)
  beta_mat_rf = normalize_rows(beta_mat)

  ## number of wells with each alpha or beta
  wa_all = Matrix::colSums(alpha_mat != 0)
  wb_all = Matrix::colSums(beta_mat != 0)
  ## total read counts of each alpha or beta
  rc_alpha = Matrix::colSums(alpha_mat)
  rc_beta = Matrix::colSums(beta_mat)
  ## total read frequency of each alpha or beta
  rf_alpha = Matrix::colSums(alpha_mat_rf)
  rf_beta = Matrix::colSums(beta_mat_rf)

  if(normalize_counts) {
    alpha_mat = alpha_mat_rf
    beta_mat = beta_mat_rf
  }

  is_lazy = !"data.frame" %in% class(well_data$alpha_meta)
  if(is_lazy) {
    if(chain == "alpha") {
      idx = find_row_index(well_data$alpha_meta_file, "targetSequences", nuc_seq)[1]
      if(is.na(idx)) stop("Nucleotide sequence not found")
      #wa1 = wa_all[idx]
      alpha_meta = extract_rows_by_index_parquet(well_data$alpha_meta_file, idx)
    } else if(chain == "beta") {
      idx = find_row_index(well_data$beta_meta_file, "targetSequences", nuc_seq)[1]
      if(is.na(idx)) stop("Nucleotide sequence not found")
      #wb1 = wb_all[idx]
      beta_meta = extract_rows_by_index_parquet(well_data$beta_meta_file, idx)
    }
  } else {
    if(chain == "alpha") {
      idx = which(well_data$alpha_meta$targetSequences == nuc_seq)[1]
      #wa1 = wa_all[idx]
      alpha_meta = well_data$alpha_meta[idx,]
    } else if(chain == "beta") {
      idx = which(well_data$beta_meta$targetSequences == nuc_seq)[1]
      if(is.na(idx)) stop("Nucleotide sequence not found")
      #wb1 = wb_all[idx]
      beta_meta = well_data$beta_meta[idx,]
    }
  }

  if(chain == "alpha") {
    vec = extract_col_dense(alpha_mat, idx)
    cor_df = sparse_col_cor(beta_mat, vec)
    cor_df$wij = sparse_overlap(beta_mat, vec)
    cor_df$wa = wa_all[idx]
    cor_df$wb = wb_all
    cor_df$readCount_alpha = rc_alpha[idx]
    cor_df$readCount_beta = rc_beta
    cor_df$readFraction_alpha = rf_alpha[idx]
    cor_df$readFraction_beta = rf_beta
  } else if(chain == "beta"){
    vec = extract_col_dense(beta_mat, idx)
    cor_df = sparse_col_cor(alpha_mat, vec)
    cor_df$wij = sparse_overlap(alpha_mat, vec)
    cor_df$wa = wa_all
    cor_df$wb = wb_all[idx]
    cor_df$readCount_alpha = rc_alpha
    cor_df$readCount_beta = rc_beta[idx]
    cor_df$readFraction_alpha = rf_alpha
    cor_df$readFraction_beta = rf_beta[idx]
  }

  cor_df$wi = cor_df$wa - cor_df$wij
  cor_df$wj = cor_df$wb - cor_df$wij

  cor_df$p_adj = p.adjust(cor_df$p, method = "BH")

  keep_indices = which(cor_df$wij > wij_cutoff*max(cor_df$wij))
  cor_df_sub = cor_df[keep_indices,]

  cor_df_sub$fisher_pval = sapply(1:nrow(cor_df_sub), function(i) {
    m11 = cor_df_sub$wij[i]
    m12 = cor_df_sub$wi[i]
    m21 = cor_df_sub$wj[i]
    m22 = length(vec) - m11 - m12 - m21
    mat = matrix(c(m11, m12, m21, m22), nrow = 2)
    tmp = fisher.test(mat)
    return(tmp$p.value)
  })

  if(is_lazy) {
    if(chain == "alpha") other_meta_sub = extract_rows_by_index_parquet(well_data$beta_meta_file, row_indices = keep_indices)
    if(chain == "beta") other_meta_sub = extract_rows_by_index_parquet(well_data$alpha_meta_file, row_indices = keep_indices)
  } else {
    if(chain == "alpha") other_meta_sub = well_data$beta_meta[keep_indices,]
    if(chain == "beta") other_meta_sub = well_data$alpha_meta[keep_indices,]
  }

  #print(other_meta_sub)

  cor_df_sub = bind_cols(other_meta_sub, cor_df_sub) %>% as_tibble() %>%
    select(chain, row_id, aaSeqCDR3, r, t, p, p_adj, wa, wb, wij, n, targetSequences, v,j, everything()) %>%
    arrange(desc(r))


  cor_df_top_n = cor_df_sub[1:n_plot,]

  # meta_gg = extract_rows_by_index_parquet(well_data$beta_meta_file, cor_df_top_n$row_id) %>%
  #   mutate(chain_name = paste("beta", row_id))
  meta_gg = cor_df_top_n %>% mutate(chain_name = paste(chain, row_id))
  if(chain == "alpha") mat_gg = as.matrix(beta_mat[,cor_df_top_n$row_id])
  if(chain == "beta") mat_gg = as.matrix(alpha_mat[,cor_df_top_n$row_id])

  colnames(mat_gg) = meta_gg$targetSequences
  rownames(mat_gg) = well_data$col_meta$well
  gg_df = lapply(1:n_plot, function(i) {
    df_tmp = tibble(readFraction1 = mat_gg[,i],
                    readFraction2 = vec,
                    well = well_data$col_meta$well,
                    chain_name = meta_gg$chain_name[i]
    )
  }) %>% bind_rows()
  gg_df$chain_name = factor(gg_df$chain_name, levels = meta_gg$chain_name)
  gg_df$readFraction1[gg_df$readFraction1 == 0] = NA
  #print(cor_df_top_n)
  gg = ggplot(gg_df, aes(x=readFraction1, y=readFraction2, label = well)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~chain_name, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  #print(gg)
  gg2 = ggplot(cor_df_sub, aes(x=r, y=-log10(p), aa = aaSeqCDR3, v = v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  gg3 = ggplot(cor_df_sub, aes(x=r, y=-log10(p_adj), aa = aaSeqCDR3, v= v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  if(interactive) {
    gg = plotly::ggplotly(gg)
    gg2 = plotly::ggplotly(gg2)
    gg3 = plotly::ggplotly(gg3)
  }

  return(list(cor_df = cor_df, cor_df_sub = cor_df_sub, cor_df_top_n = cor_df_top_n, gg = gg, gg2 = gg2, gg3 = gg3))
}
