seq_error_detect_test = function(well_data, chain, normalize_counts = TRUE, dist_cutoff = 4, p_cutoff = 0.05) {
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
  rf_alpha = rc_alpha/sum(rc_alpha)
  rf_beta = rc_beta/sum(rc_beta)

  if(normalize_counts) {
    alpha_mat = alpha_mat_rf
    beta_mat = beta_mat_rf
  }

  alpha_mat_log = alpha_mat_rf
  alpha_mat_log@x = log10(alpha_mat_log@x)
  beta_mat_log = beta_mat_rf
  beta_mat_log@x = log10(beta_mat_log@x)

  mat = get(paste(chain, "mat", sep = "_"))

  #is_lazy = !"data.frame" %in% class(well_data$alpha_meta)
  vj = paste( pull(well_data[[paste(chain, "meta", sep = "_")]], v), pull(well_data[[paste(chain, "meta", sep = "_")]], j) )
  vj_tbl = table(vj) %>% sort() %>% rev()
  vj_tbl2 = vj_tbl[which(vj_tbl > 1)]
  ## for each VJ combo, take only those clones and get pairwise correlation + edit distance
  for(i in 1:length(vj_tbl2)) {
    idx = which(vj == names(vj_tbl2)[i])
    #beta_meta_sub = well_data$beta_meta[idx,]
    meta_sub = extract_rows_by_index_parquet(well_data[[paste(chain, "_meta_file", sep = "")]], idx)
    mat_sub = mat[,idx]
    cor_mat = sparse_col_cor_all(mat_sub, mat_sub)
    ## for each clone from most abundant to least, get other clones with p<=p_cutoff and dist <= dist_cutoff (r>0) and discard them
    ## report how many discarded
    for(j in 1:ncol(mat_sub)) {
      meta_sub$cor = cor_list$correlation[,j]
      meta_sub$pval = cor_list$p_value[,j]
      meta_sub$tstat = cor_list$t_statistic[,j]

    }

  }

  ### testing
  #ll = sparse_col_cor_all(alpha_mat[,1:2], beta_mat)
  #cor_df = sparse_col_cor(beta_mat[,(i+1):ncol(beta_mat)], beta_mat[,i])
  cor_df = sparse_col_cor(beta_mat, beta_mat[,i])
  cor_df$readCount = rc_beta
  cor_df$readFraction = rf_beta
  cor_df$vj = paste(well_data$beta_meta$v, well_data$beta_meta$j)
  cor_df$p_adj = p.adjust(cor_df$p, method = "BH")
  test = bind_cols(well_data$beta_meta,cor_df) %>%
    arrange(desc(r)) %>% as_tibble()
  test$dist_hamming = stringdist::stringdist(well_data$beta_meta$targetSequences[1], test$targetSequences, method = "hamming", nthread = 4)
  test$dist_edit = stringdist::stringdist(well_data$beta_meta$targetSequences[1], test$targetSequences, method = "lv", nthread = data.table::getDTthreads())

  test2 = test[-1,]
  test_cluster = test2 %>% filter(p <= p_cutoff, dist_edit <= dist_cutoff, r > 0)
  test_keep = test2 %>% filter(p > p_cutoff | dist_edit > dist_cutoff | r <= 0)

  hist(test$p)

  hist(test_keep %>% filter(dist_edit < 20) %>% extract2("dist_edit"))

}
