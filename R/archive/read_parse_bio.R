

read_parse_bio = function(path) {
  dt_parse<-fread(path)
  dt_parse_clean<-get_clonotypes_parse_bio(dt_parse) %>% as.data.table()
  ## dplyr version
  # df_complete = df_test %>%
  #   arrange(desc(n_cells)) %>%
  #   filter(!duplicated(paste0(beta_nuc,"_",alpha_nuc))) %>%
  #   filter(!is.na(beta_nuc)&!is.na(alpha_nuc))
  ## data.table version
  dt_parse_complete<-dt_parse_clean[order(-n_cells),][!duplicated(paste0(beta_nuc,"_",alpha_nuc)),][!is.na(beta_nuc)&!is.na(alpha_nuc),]
  #dt_parse_complete[,alpha_beta:=paste0(alpha_nuc,"_",beta_nuc),]
  #dt_parse_complete[,beta_nuc:=beta_nuc,]
  #dt_parse_complete[,alpha_nuc:=alpha_nuc,]
  #' note: n_cells contains extremely repeated values, so there are often many
  #' rows tied with highest "n_cells" count for each unique alpha_beta.
  #' Selecting the row with the highest "n_cells" count for each makes this step
  #' non-repeatable bc of different ordering of ties.
  #' This causes occasional mismatches in barcode, va, cdr3a2, va2, ja2
  list(complete=dt_parse_complete,clean=dt_parse_clean,raw=dt_parse)
}


read_parse_bio_tcr_annotation_airr = function(path) {

}

#df = fread("~/git/newell_benchmarking/data/BM02_Parse_v2_F5_25k/TCR/COMBINED/B1-PA-F5/TCR_filtered/tcr_annotation_airr.tsv")
get_clonotypes_parse_bio = function(df) {
  ## make df with top two alphas for each cell barcode (most umis)
  df_alpha = df %>% filter(locus == "TRA")  %>% group_by(cell_barcode) %>% arrange(desc(transcript_count)) %>%
    mutate(alpha_id = row_number()) %>% slice_head(n=2) %>%
    select(cell_barcode, cdr3_aa, cdr3, v_call, j_call, read_count, transcript_count, alpha_id) %>%
    ungroup()
  ## make df with first alpha for each cell barcode
  df_alpha1 = df_alpha %>% filter(alpha_id == 1) %>% ## first alpha chain
    dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3, va = v_call, ja = j_call,
                  reads_alpha = read_count, umis_alpha = transcript_count) %>%
    select(-alpha_id) %>%
    ungroup()
  ## make df with second alpha for each cell barcode
  df_alpha2 = df_alpha %>% filter(alpha_id == 2) %>% ## second alpha chain
    dplyr::rename(cdr3a2 = cdr3_aa, cdr3a2_nt = cdr3, va2 = v_call, ja2 = j_call,
                  reads_alpha2 = read_count, umis_alpha2 = transcript_count) %>%
    select(-alpha_id) %>%
    ungroup()
  ## make df with top beta for each cell barcode (most umis)
  df_beta = df %>% filter(locus == "TRB")  %>% group_by(cell_barcode) %>% arrange(desc(transcript_count)) %>%
    mutate(beta_id = row_number()) %>% slice_head(n=1) %>%
    select(cell_barcode, cdr3_aa, cdr3, v_call, j_call, read_count, transcript_count) %>%
    dplyr::rename(cdr3b = cdr3_aa, cdr3b_nt = cdr3, vb = v_call, jb = j_call,
                  reads_beta = read_count, umis_beta = transcript_count) %>%
    ungroup()
  ## join dfs to make df with top beta and top two alphas for each cell barcode
  df_all = df_beta %>% full_join(df_alpha1, by=c("cell_barcode")) %>% full_join(df_alpha2, by=c("cell_barcode")) %>%
    ungroup() %>%
    mutate(alpha_beta = paste(cdr3a_nt, cdr3b_nt, sep = "_")) %>%
    add_count(alpha_beta, name = "n_cells") %>%
    mutate(has_beta = !is.na(cdr3b_nt), has_alpha = !is.na(cdr3a_nt), has_alpha2 = !is.na(cdr3a2_nt)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    arrange(desc(n_cells))  %>%
    dplyr::rename(beta_nuc = cdr3b_nt, alpha_nuc = cdr3a_nt, alpha2_nuc = cdr3a2_nt)
  return(df_all)
}
