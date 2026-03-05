#' Calculate and plot read fraction correlation across wells (T-SHELL)
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function takes in an object with well-level read counts for a sample,
#' loaded by \code{\link{load_well_counts_binary}()} along with a CDR3 sequence
#' of interest for TCRalpha or TCRbeta. It calculates the correlation of the input
#' sequence with all of its potential partner chains and returns the top N chains
#' with the highest correlation/T-SHELL value.
#'
#' @param well_data well-level read count data loaded by \code{\link{load_well_counts_binary}()}
#' @param nuc_seq the CDR3 nucleotide sequence of the chain of interest
#' @param chain the chain of the input nucleotide sequence, either "alpha" or "beta"
#' @param n_plot the number of scatter plots to return (one for each of top N
#' potential partners, default is 9)
#' @param plot whether to plot the data or only return data frames (default is TRUE)
#' @param plot_log_scale whether to plot log(readFraction) instead of readFraction (default is FALSE)
#' @param interactive whether to return interactive plots (default is FALSE)
#' @param normalize_counts whether to normalize the counts by the total count for each well (default is TRUE)
#' @param loss_frac_cutoff the cutoff to use for "loss fraction" of potential partners for plots. Must be between 0 and 1. 1 or higher will not filter any clones, lower values will apply more strict filtering. (default is 0.5).
#' @param calc_fisher_pval_all whether to calculate the Fisher exact test p-value
#' of well overlap between the input chain and all of the filtered clones.
#' Default is FALSE because this is computationally demanding.
#'
#' @return
#' A list of output with the following elements:
#' - data (list)
#'  - df_top_n - a data frame with results for the top N potential partner chains
#'  - df_all - a data frame with results for all potential partners passing the \code{loss_frac_cutoff}
#'  - input_meta - a one-row data frame with the metadata for the input chain
#' - plots (list)
#'  - scatter - a scatter plot of read fractions for the input chain vs. the
#'  top N potential partner chains
#'  - r_vs_p - a point plot of correlation (r) vs. -log10(raw p-value)
#'  - r_vs_p_adj - a point plot of correlation (r) vs. -log10(FDR adjusted p-value)
#'  - rank_vs_p_adj - a Manhattan plot of chain rank vs. -log10(FDR adjusted p-value)
#' - call - a list with all arguments given to the function
#'
#' @family well
#'
#' @export
#'

plot_tshell = function(
    well_data,
    nuc_seq,
    chain,
    n_plot = 9,
    plot = TRUE,
    plot_log_scale = FALSE,
    interactive = FALSE,
    normalize_counts = TRUE,
    loss_frac_cutoff = 0.5, ## change to 1 or Inf to ignore filter
    calc_fisher_pval_all = FALSE
) {
  call = as.list(environment())
  call$well_data = NULL
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

  alpha_mat_log = alpha_mat_rf
  alpha_mat_log@x = log10(alpha_mat_log@x)
  beta_mat_log = beta_mat_rf
  beta_mat_log@x = log10(beta_mat_log@x)


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
    vec_meta = alpha_meta
    cor_df = sparse_col_cor(beta_mat, vec)
    cor_df_log_log = sparse_col_cor(beta_mat_log, extract_col_dense(alpha_mat_log, idx))
    cor_df$r_log_log = cor_df_log_log$r
    cor_df$p_log_log = cor_df_log_log$p
    cor_df$wij = sparse_overlap(beta_mat, vec)
    cor_df$wa = wa_all[idx]
    cor_df$wb = wb_all
    cor_df$readCount_alpha = rc_alpha[idx]
    cor_df$readCount_beta = rc_beta
    cor_df$readFraction_alpha = rf_alpha[idx]
    cor_df$readFraction_beta = rf_beta
  } else if(chain == "beta"){
    vec = extract_col_dense(beta_mat, idx)
    vec_meta = beta_meta
    cor_df = sparse_col_cor(alpha_mat, vec)
    cor_df_log_log = sparse_col_cor(alpha_mat_log, extract_col_dense(beta_mat_log, idx))
    cor_df$r_log_log = cor_df_log_log$r
    cor_df$p_log_log = cor_df_log_log$p
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
  cor_df$p_adj_old = cor_df$p/(sort(cor_df$p)[3])

  cor_df$loss_a_frac = with(cor_df, (wb-wij)/(wij+(wb-wij)+(wa-wij)))
  cor_df$loss_b_frac = with(cor_df, (wa-wij)/(wij+(wb-wij)+(wa-wij)))
  cor_df$loss_frac_sum = with(cor_df, loss_a_frac+loss_b_frac)
  #keep_indices = which(cor_df$wij > wij_cutoff*max(cor_df$wij))
  keep_indices = which(cor_df$loss_frac_sum <= loss_frac_cutoff)

  cor_df_sub = cor_df[keep_indices,]

  if(calc_fisher_pval_all) cor_df_sub$fisher_pval = calc_fisher_pval_df(cor_df_sub, n= length(vec))

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
    arrange(desc(t)) %>% mutate(rank = row_number())


  nrow_final = min(n_plot, nrow(cor_df_sub))
  cor_df_top_n = cor_df_sub[1:nrow_final,] %>%
    mutate(chain_name = paste(paste("rank ", rank, ":", sep = ""), chain, row_id, "\n", aaSeqCDR3)) %>%
    mutate(chain_name = factor(chain_name, levels = chain_name))

  if(!calc_fisher_pval_all) cor_df_top_n$fisher_pval = calc_fisher_pval_df(cor_df_top_n, n =length(vec))

  # meta_gg = extract_rows_by_index_parquet(well_data$beta_meta_file, cor_df_top_n$row_id) %>%
  #   mutate(chain_name = paste("beta", row_id))
  meta_gg = cor_df_top_n
  if(chain == "alpha") mat_gg = as.matrix(beta_mat[,cor_df_top_n$row_id])
  if(chain == "beta") mat_gg = as.matrix(alpha_mat[,cor_df_top_n$row_id])

  colnames(mat_gg) = meta_gg$targetSequences
  rownames(mat_gg) = well_data$well_meta$well
  gg_df = lapply(1:nrow_final, function(i) {
    df_tmp = tibble(readFraction1 = mat_gg[,i],
                    readFraction2 = vec,
                    well = well_data$well_meta$well,
                    chain_name = meta_gg$chain_name[i]
    )
  }) %>% bind_rows()
  #gg_df$chain_name = factor(gg_df$chain_name, levels = meta_gg$chain_name)
  gg_df$readFraction1[gg_df$readFraction1 == 0] = NA
  #print(cor_df_top_n)
  label_df = cor_df_top_n %>%
    mutate(
      #label = sprintf("r = %.2f\np_adj = ", r, scales::scientific(p_adj)),
      label   = paste0("r = ", round(r, 2),"\np_adj = ", scales::scientific(p_adj))#,
      #chain_name = paste(chain, row_id)
      ) #%>%
    #mutate(chain_name = factor(chain_name, levels = meta_gg$chain_name))
  if(plot) {
    gg = ggplot(gg_df, aes(x=readFraction1, y=readFraction2, label = well)) +
      geom_point(alpha = 0.3) +
      geom_text(data = label_df, aes(x = 0, y = Inf, label = label),
                hjust = -0.1, vjust = 1.5, size = 3, inherit.aes = FALSE) +
      facet_wrap(~chain_name, scales = "free") +
      #scale_x_continuous(labels = scales::label_scientific()) +
      #scale_y_continuous(labels = scales::label_scientific()) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if(plot_log_scale) {
      gg = gg + scale_x_log10() + scale_y_log10()
    }
    #print(gg)
    gg2 = ggplot(cor_df_sub, aes(x=r, y=-log10(p), aa = aaSeqCDR3, v = v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
    gg3 = ggplot(cor_df_sub, aes(x=r, y=-log10(p_adj), aa = aaSeqCDR3, v= v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
    gg4 = ggplot(cor_df_sub, aes(x = rank, y = -log10(p_adj))) +
      geom_segment(aes(x = rank, xend = rank, y = 0, yend = -log10(p_adj))) +          # color by significance
      geom_hline(yintercept = -log10(0.05),          # significance threshold line
                 linetype = "dashed", color = "red") +
      geom_point(data = filter(cor_df_sub, rank == 1), size = 4, color = "darkgreen") +
      geom_text(data = label_df[1,], aes(x = 1, y = Inf, label = label),
                hjust = -0.5, vjust = 1.5, size = 3, inherit.aes = FALSE) +
      labs(x = "Rank", y = expression(-log[10](p))) +
      theme_minimal()
    if(interactive) {
      gg = plotly::ggplotly(gg)
      gg2 = plotly::ggplotly(gg2)
      gg3 = plotly::ggplotly(gg3)
      gg4 = plotly::ggplotly(gg4)
    }
  } else {
    gg = NULL
    gg2 = NULL
    gg3 = NULL
    gg4 = NULL
  }


  return(
    list(data = list(input_meta = vec_meta, df_all = cor_df_sub, df_top_n = cor_df_top_n),
         plots = list(scatter = gg, rank_vs_p_adj = gg4, r_vs_p = gg2, r_vs_p_adj = gg3),
         call = call
         )
    )
}

