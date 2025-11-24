# Imputes missing alpha and beta chains where possible for single-cell TIRTLseq data
#
# @description
# `r lifecycle::badge('experimental')`
#
# @param df a data frame with single-cell TIRTLseq data
# @param verbose print a before and after summary of clonotypes (default TRUE)
# @param keep_all_columns keep intermediate columns created by the function (default FALSE)
#
# @family single-cell
#
clean_scTIRTLseq = function(df, verbose = TRUE, keep_all_columns = FALSE) {
  df = df %>%
    mutate(id = 1:nrow(df))
  df_orig = df
  df_copy = df %>%
    mutate(across(c(aaSeqCDR3u_beta, bestV_beta, bestJ_beta, nSeqCDR3_beta,
                    aaSeqCDR3u_alpha, bestV_alpha, bestJ_alpha, nSeqCDR3_alpha,
                    aaSeqCDR3u_alpha_second, bestV_alpha_second, bestJ_alpha_second, nSeqCDR3_alpha_second
    ), ~tidyr::replace_na(., ""), .names = "noNA_{.col}")) %>%
    mutate(beta_full = paste(noNA_aaSeqCDR3u_beta, noNA_bestV_beta, noNA_bestJ_beta, noNA_nSeqCDR3_beta, sep = "|")) %>%
    mutate(alpha_full = paste(noNA_aaSeqCDR3u_alpha, noNA_bestV_alpha, noNA_bestJ_alpha, noNA_nSeqCDR3_alpha, sep = "|")) %>%
    mutate(alpha_full2 = paste(noNA_aaSeqCDR3u_alpha_second, noNA_bestV_alpha_second, noNA_bestJ_alpha_second, noNA_nSeqCDR3_alpha_second, sep = "|")) %>%
    mutate(
      across(c(alpha_full, alpha_full2, beta_full,
               aaSeqCDR3u_beta, bestV_beta, bestJ_beta, nSeqCDR3_beta,
               aaSeqCDR3u_alpha, bestV_alpha, bestJ_alpha, nSeqCDR3_alpha,
               aaSeqCDR3u_alpha_second, bestV_alpha_second, bestJ_alpha_second, nSeqCDR3_alpha_second
      ), identity, .names = "{.col}_orig")
    )
  ## fill in missing alphas
  group_var = "beta_full"
  df_beta = df_copy %>%
    filter(alpha_full != "|||")
  betas = table(df_copy[[group_var]]) %>% sort() %>% rev() %>% names()
  betas = betas[betas != "|||"]
  df_beta2 = lapply(betas, function(x) {
    df_tmp = df_beta[df_beta[[group_var]] == x,]
    tab = rev(sort(table(c(df_tmp$alpha_full, df_tmp$alpha_full2))))
    tab = tab[names(tab) != "|||"]
    df_ret = tibble(beta = x, alpha1 = names(tab)[1], alpha2 = names(tab)[2], n_alpha1 = tab[1], n_alpha2 = tab[2], n_beta = dim(df_tmp)[1], alpha_3 = names(tab)[3], n_alpha3 = tab[3])
  }) %>% bind_rows() %>%
    arrange(desc(n_beta)) %>%
    filter(!is.na(alpha1))

  ## fill in missing alphas
  group_var = "alpha_full"
  group_var2 = "alpha_full2"
  df_alpha = df_copy %>%
    filter(beta_full != "|||")
  alphas = table( c(df_alpha[[group_var]], df_alpha[[group_var2]]) ) %>% sort() %>% rev() %>% names()
  alphas = alphas[alphas != "|||"]
  df_alpha2 = lapply(alphas, function(x) {
    df_tmp1 = df_alpha[df_alpha[[group_var]] == x,]
    df_tmp2 = df_alpha[df_alpha[[group_var2]] == x,]
    df_tmp = bind_rows(df_tmp1, df_tmp2)
    tab = rev(sort(table(df_tmp$beta_full)))
    tab = tab[names(tab) != "|||"]
    df_ret = tibble(alpha = x, beta1 = names(tab)[1], beta2 = names(tab)[2], n_beta1 = tab[1], n_beta2 = tab[2], n_alpha = dim(df_tmp)[1], beta_3 = names(tab)[3], n_beta3 = tab[3])
  }) %>% bind_rows() %>%
    arrange(desc(n_alpha)) %>%
    filter(!is.na(beta1))


  df_fill_beta = lapply(1:length(df_beta2$beta), function(i) {
    #df_tmp_orig = df_copy[df_copy$beta_full == df_beta2$beta[i],]
    df_tmp = df_copy[df_copy$beta_full == df_beta2$beta[i],]
    alpha1_full = df_beta2$alpha1[i]
    alpha2_full = df_beta2$alpha2[i]
    df_tmp$change = ""
    for(j in 1:nrow(df_tmp)) {
      if(is.na(alpha2_full)) {
        if(df_tmp$alpha_full[j] == alpha1_full) {
          df_tmp$change[j] = "no change -- alpha observed"
        } else {
          alpha1_vec = strsplit(alpha1_full, "\\|")[[1]]
          df_tmp$alpha_full[j] = alpha1_full
          df_tmp$aaSeqCDR3u_alpha[j] = alpha1_vec[1]
          df_tmp$bestV_alpha[j] = alpha1_vec[2]
          df_tmp$bestJ_alpha[j] = alpha1_vec[3]
          df_tmp$nSeqCDR3_alpha[j] = alpha1_vec[4]
          df_tmp$change[j] = "added alpha"
        }
      } else {
        if(df_tmp$alpha_full[j] == alpha1_full && df_tmp$alpha_full2[j] == alpha2_full) {
          df_tmp$change = "no change -- both alphas observed"
        } else if(df_tmp$alpha_full[j] == alpha2_full && df_tmp$alpha_full2[j] == alpha1_full) {
          df_tmp$change = "switched alpha1 and alpha2"
        } else if(df_tmp$alpha_full[j] == alpha1_full && df_tmp$alpha_full2[j] == "|||") {
          df_tmp$change = "added alpha2"
        } else if(df_tmp$alpha_full[j] == alpha2_full && df_tmp$alpha_full2[j] == "|||"){
          df_tmp$change = "added alpha1"
        } else if(df_tmp$alpha_full[j] == "|||") {
          df_tmp$change = "added both alphas"
        } else {
          df_tmp$change = "other"
        }
        alpha1_vec = strsplit(alpha1_full, "\\|")[[1]]
        df_tmp$alpha_full[j] = alpha1_full
        df_tmp$aaSeqCDR3u_alpha[j] = alpha1_vec[1]
        df_tmp$bestV_alpha[j] = alpha1_vec[2]
        df_tmp$bestJ_alpha[j] = alpha1_vec[3]
        df_tmp$nSeqCDR3_alpha[j] = alpha1_vec[4]
        alpha2_vec = strsplit(alpha2_full, "\\|")[[1]]
        df_tmp$alpha_full2[j] = alpha2_full
        df_tmp$aaSeqCDR3u_alpha_second[j] = alpha2_vec[1]
        df_tmp$bestV_alpha_second[j] = alpha2_vec[2]
        df_tmp$bestJ_alpha_second[j] = alpha2_vec[3]
        df_tmp$nSeqCDR3_alpha_second[j] = alpha2_vec[4]
      }
    }
    return(df_tmp)
  }) %>% bind_rows()

  df_fill_alpha = lapply(1:length(df_alpha2$alpha), function(i) {
    df_tmp = df_copy[(df_copy$alpha_full == df_alpha2$alpha[i]) | (df_copy$alpha_full2 == df_alpha2$alpha[i]),]
    beta1_full = df_alpha2$beta1[i]
    beta2_full = df_alpha2$beta2[i]
    df_tmp$change = ""
    for(j in 1:nrow(df_tmp)) {
      if(is.na(beta2_full)) {
        if(df_tmp$beta_full[j] == beta1_full) {
          df_tmp$change[j] = "no change -- beta observed"
        } else {
          df_tmp$change[j] = "added beta"
          beta1_vec = strsplit(beta1_full, "\\|")[[1]]
          df_tmp$beta_full[j] = beta1_full
          df_tmp$aaSeqCDR3u_beta[j] = beta1_vec[1]
          df_tmp$bestV_beta[j] = beta1_vec[2]
          df_tmp$bestJ_beta[j] = beta1_vec[3]
          df_tmp$nSeqCDR3_beta[j] = beta1_vec[4]
        }
      } else {
        if(df_tmp$beta_full[j] == beta1_full) {
          df_tmp$change[j] = "no change -- possible second beta"
        } else if(df_tmp$beta_full[j] == beta2_full) {
          df_tmp$change[j] = "no change -- possible second beta"
        } else {
          df_tmp$change[j] = "added beta -- possible second beta"
          beta1_vec = strsplit(beta1_full, "\\|")[[1]]
          df_tmp$beta_full[j] = beta1_full
          df_tmp$aaSeqCDR3u_beta[j] = beta1_vec[1]
          df_tmp$bestV_beta[j] = beta1_vec[2]
          df_tmp$bestJ_beta[j] = beta1_vec[3]
          df_tmp$nSeqCDR3_beta[j] = beta1_vec[4]
        }
      }
    }
    return(df_tmp)
  }) %>% bind_rows()

  df_fill_beta2 = df_fill_beta %>% filter(!grepl("no change", change))
  df_fill_alpha2 = df_fill_alpha %>% filter(!grepl("no change", change))

  df_imputed = bind_rows(df_fill_alpha2, df_fill_beta2)
  df_return = bind_rows(df_copy %>% filter(!id %in% df_imputed$id), df_imputed) %>%
    arrange(id) %>%
    mutate_at(vars(aaSeqCDR3u_beta, bestV_beta, bestJ_beta, nSeqCDR3_beta,
                   aaSeqCDR3u_alpha, bestV_alpha, bestJ_alpha, nSeqCDR3_alpha), ~replace(., .=="" , NA)) %>%
    mutate_at(vars(alpha_full, alpha_full2, beta_full), ~replace(., .=="|||" , NA))
  if(!keep_all_columns) {
    df_return = df_return %>%
      select(!starts_with("noNA_")) %>%
      select(!ends_with("_orig"))
  } else {
    df_return = df_return %>%
      select(c(!starts_with("noNA_"), !ends_with("_orig")), ends_with("_orig"), starts_with("noNA_"))
  }

  if(verbose) {
    message("Before cleaning/imputation:")
    summarize_scTIRTLseq(df_orig)

    message("After cleaning/imputation:")
    summarize_scTIRTLseq(df_return)
  }

  return(list(df_new = df_return, df_orig = df_orig, clones_by_alpha = df_alpha2, clones_by_beta = df_beta2))
}
