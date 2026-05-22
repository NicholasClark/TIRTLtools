#df = fread("~/git/newell_benchmarking/data/BM02_GEMX5P_F5_25k/multi/count/")
#df = fread("~/git/newell_benchmarking/data/GSM/GSM8951948_GEMX5P_F5_filtered_contig_annotations.csv.gz")
#df = fread("~/git/newell_benchmarking/data/GSM/GSM8951952_NextGEM5P_F5_filtered_contig_annotations.csv.gz")

read_10x = function(path) {
  dt_10x<-fread(path)
  dt_10x_clean<-get_clonotypes_10x(dt_10x) %>% as.data.table()
  ## dplyr version
  # df_complete = df_test %>%
  #   arrange(desc(n_cells)) %>%
  #   filter(!duplicated(paste0(beta_nuc,"_",alpha_nuc))) %>%
  #   filter(!is.na(beta_nuc)&!is.na(alpha_nuc))
  ## data.table version
  dt_10x_complete<-dt_10x_clean[order(-n_cells),][!duplicated(paste0(beta_nuc,"_",alpha_nuc)),][!is.na(beta_nuc)&!is.na(alpha_nuc),]

  list(complete=dt_10x_complete,clean=dt_10x_clean,raw=dt_10x)
}


get_clonotypes_10x = function(df, allow_second_alpha = TRUE, separate_row_for_second_alpha = TRUE) {
  ## make df with top two alphas for each cell barcode (most umis)
  df_alpha_all = df %>% filter(chain == "TRA")  %>% group_by(barcode) %>% arrange(desc(umis), desc(reads), .by_group = T) %>%
    mutate(alpha_id = row_number(), n_alphas = n()) %>%
    select(barcode, cdr3, cdr3_nt, v_gene, j_gene, reads, umis, alpha_id, n_alphas) %>%
    ungroup()
  ## make df with first alpha for each cell barcode
  df_alpha1 = df_alpha_all %>% filter(alpha_id == 1) %>% ## first alpha chain
    dplyr::rename(cdr3a = cdr3, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis) %>%
    select(-alpha_id) %>%
    ungroup()

  df_alpha_multi = df_alpha_all %>% filter(alpha_id > 2) %>% ## second alpha chain
    dplyr::rename(cdr3a2 = cdr3, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis) %>%
    select(-alpha_id) %>%
    ungroup()
  ## make df with top beta for each cell barcode (most umis)
  df_beta_all = df %>% filter(chain == "TRB")  %>% group_by(barcode) %>% arrange(desc(umis), .by_group = T) %>%
    mutate(beta_id = row_number(), n_betas = n()) %>%
    select(barcode, cdr3, cdr3_nt, v_gene, j_gene, reads, umis, beta_id, n_betas) %>%
    dplyr::rename(cdr3b = cdr3, cdr3b_nt = cdr3_nt, vb = v_gene, jb = j_gene, reads_beta = reads, umis_beta = umis) %>%
    ungroup()
  df_beta = df_beta_all %>% filter(beta_id == 1)
  df_beta_multi = df_beta_all %>% filter(beta_id > 1)

  df_beta_test = df_beta_all %>% filter(n_betas > 1) %>% arrange(barcode)
  df_alpha_test = df_alpha_all %>% filter(n_alphas > 1) %>% arrange(barcode)

  ## join dfs to make df with top beta and top two alphas for each cell barcode

  ## make df with second alpha for each cell barcode
  if(allow_second_alpha) {
    df_alpha2 = df_alpha_all %>% filter(alpha_id == 2) %>% ## second alpha chain
      select(-alpha_id) %>%
      ungroup()
    if(!separate_row_for_second_alpha) {
      df_alpha2 = df_alpha2 %>% dplyr::rename(cdr3a2 = cdr3, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis)
      df_all1 = df_beta %>%
        full_join(df_alpha1, by=c("barcode")) %>%
        full_join(df_alpha2, by=c("barcode")) %>%
        ungroup()
    } else {
      df_alpha2 = df_alpha2 %>% dplyr::rename(cdr3a = cdr3, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_alpha12 = bind_rows(df_alpha1, df_alpha2)
      df_all1 = df_beta %>%
        full_join(df_alpha12, by=c("barcode")) %>%
        ungroup()
    }
  } else {
    df_all1 = df_beta %>%
      full_join(df_alpha1, by=c("barcode")) %>%
      ungroup()
  }

  df_all = df_all1 %>%
    mutate(alpha_beta = paste(cdr3a_nt, cdr3b_nt, sep = "_")) %>%
    add_count(alpha_beta, name = "n_cells") %>%
    mutate(has_beta = !is.na(cdr3b_nt), has_alpha = !is.na(cdr3a_nt)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    arrange(desc(n_cells)) %>%
    dplyr::rename(beta_nuc = cdr3b_nt, alpha_nuc = cdr3a_nt)
  if(allow_second_alpha && !separate_row_for_second_alpha) {
    df_all = df_all %>% mutate(has_alpha2 = !is.na(cdr3a2_nt)) %>% dplyr::rename(alpha2_nuc = cdr3a2_nt)
  }

  return(df_all)
}
