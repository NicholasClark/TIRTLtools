read_external = function(path, format = c("auto","10X", "ParseBio"), multi = TRUE, separate_rows = TRUE) {
  format = format[1]
  df_orig = fread(path)
  if(format == "auto") format = infer_format(df_orig)
  df = rename_columns(df_orig, format=format, rename_df = get_names_df_single_chain(), verbose = TRUE)
  df_pairs = convert_chains_to_pairs(df, multi = multi, separate_rows = separate_rows) %>% as.data.table()
  df_complete <- df_pairs[order(-n_cells),][!duplicated(paste0(beta_nuc,"_",alpha_nuc)),][!is.na(beta_nuc)&!is.na(alpha_nuc),]

  out = list(complete = df_complete, clean = df_pairs, raw = df_orig)
  return(out)
}

## multi -- if FALSE, only take one beta and best two alpha, else take all pairs
## separate_rows -- if TRUE, report each pair on a separate row, instead of extra columns for second alpha
convert_chains_to_pairs = function(df, multi = TRUE, separate_rows = TRUE) {
  df_alpha_all = df %>%
    filter(chain == "TRA") %>%
    group_by(barcode) %>%
    arrange(desc(umis), desc(reads), .by_group = T) %>%
    mutate(alpha_id = row_number(), n_alphas = n()) %>%
    select(barcode, cdr3_aa, cdr3_nt, v_gene, j_gene, reads, umis, alpha_id, n_alphas) %>%
    ungroup()
  ## make df with first alpha for each cell barcode
  df_alpha1 = df_alpha_all %>%
    filter(alpha_id == 1) %>% ## first alpha chain
    dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis) %>%
    select(-alpha_id) %>%
    ungroup()

  df_alpha_multi = df_alpha_all %>%
    filter(alpha_id > 2) %>% ## second alpha chain
    dplyr::rename(cdr3a2 = cdr3_aa, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis) %>%
    select(-alpha_id) %>%
    ungroup()
  ## make df with top beta for each cell barcode (most umis)
  df_beta_all = df %>%
    filter(chain == "TRB") %>%
    group_by(barcode) %>%
    arrange(desc(umis), .by_group = T) %>%
    mutate(beta_id = row_number(), n_betas = n()) %>%
    select(barcode, cdr3_aa, cdr3_nt, v_gene, j_gene, reads, umis, beta_id, n_betas) %>%
    dplyr::rename(cdr3b = cdr3_aa, cdr3b_nt = cdr3_nt, vb = v_gene, jb = j_gene, reads_beta = reads, umis_beta = umis) %>%
    ungroup()
  df_beta = df_beta_all %>% filter(beta_id == 1)
  df_beta_multi = df_beta_all %>% filter(beta_id > 1)

  df_beta_test = df_beta_all %>% filter(n_betas > 1) %>% arrange(barcode)
  df_alpha_test = df_alpha_all %>% filter(n_alphas > 1) %>% arrange(barcode)

  ## join dfs to make df with top beta and top two alphas for each cell barcode

  ## make df with second alpha for each cell barcode
  if(!separate_rows) {
    df_alpha2 = df_alpha_all %>%
      filter(alpha_id == 2) %>% ## use only second alpha chain
      select(-alpha_id) %>%
      ungroup() %>%
      dplyr::rename(cdr3a2 = cdr3, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis)
    df_all1 = df_beta %>%
      full_join(df_alpha1, by=c("barcode")) %>%
      full_join(df_alpha2, by=c("barcode")) %>%
      ungroup()
  } else { ## if separate_rows
    if(multi) {
      df_alpha2 = df_alpha_all %>%
        filter(alpha_id == 2) %>% ## use only second alpha chain
        select(-alpha_id) %>%
        ungroup() %>%
        dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_alpha12 = bind_rows(df_alpha1, df_alpha2)
      df_all1 = df_beta %>%
        full_join(df_alpha12, by=c("barcode")) %>%
        ungroup()
    } else {
      df_join = df_alpha_all %>%
        dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_all1 = df_beta_all %>%
        full_join(df_join, by=c("barcode")) %>%
        ungroup()
    }
  }
  df_all = df_all1 %>%
    mutate(alpha_beta = paste(cdr3a_nt, cdr3b_nt, sep = "_")) %>%
    add_count(alpha_beta, name = "n_cells") %>%
    mutate(has_beta = !is.na(cdr3b_nt), has_alpha = !is.na(cdr3a_nt)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    arrange(desc(n_cells)) %>%
    dplyr::rename(beta_nuc = cdr3b_nt, alpha_nuc = cdr3a_nt)
  if(!separate_rows) {
    df_all = df_all %>% mutate(has_alpha2 = !is.na(cdr3a2_nt)) %>% dplyr::rename(alpha2_nuc = cdr3a2_nt)
  }

  return(df_all)
}

infer_format = function(df, verbose = TRUE) {
  if("locus" %in% colnames(df)) {
    if(verbose) message("Reading as 'ParseBio'...")
    return("ParseBio")
  }
  if("chain" %in% colnames(df)) {
    if(verbose) message("Reading as '10X'...")
    return("10X")
  }
  if(verbose) message("Format could not be identified, trying to read as '10X'...")
  return("10X")
}
