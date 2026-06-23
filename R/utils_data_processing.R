.process_TIRTLseq_single = function(data, clean = FALSE, remove_nonfunctional = FALSE, verbose = TRUE) {
  checkmate::assert_list(data)
  checkmate::assert_logical(clean)
  checkmate::assert_logical(remove_nonfunctional)

  data = lapply(data, function(x) .identify_non_functional_seqs_single(x, remove = remove_nonfunctional))
  check1 = .check_has_paired_data(data, error = FALSE, warn = TRUE)
  check2 = .check_has_single_chain_data(data, error = FALSE, warn = TRUE)
  if((!check1) || (!check2)) {
    warning("Skipping cleaning and annotation")
    return(data)
  }

  data$alpha = data$alpha %>%
    arrange(desc(readFraction), desc(is_functional)) %>%
    mutate(rank=row_number(), id = paste("alpha", rank, sep = "_"))
  data$beta = data$beta %>%
    arrange(desc(readFraction), desc(is_functional)) %>%
    mutate(rank=row_number(), id = paste("beta", rank, sep = "_"))
  data$paired = .add_single_chain_data_simple(data, verbose = verbose)  %>%
    arrange(beta_rank, alpha_rank) %>%
    .get_receptor_component()


  data = data %>% .annotate_paired_single() ## adds "is_paired" column to single-chain data and adds "paired_alt" data frame

  if(clean) data$paired = data$paired %>% .clean_pairs_single()
  if(clean && (!is.null(data$paired_alt))) data$paired_alt = data$paired_alt %>% .clean_pairs_single()

  data$paired = data$paired %>%
    .order_columns(.get_paired_cols_ordered(), keep_extra = TRUE, verbose = verbose) %>%
    arrange(receptor_rank)
  data$paired_alt = data$paired_alt %>%
    .order_columns(.get_paired_alt_cols_ordered(), keep_extra = TRUE, verbose = verbose) %>%
    arrange(receptor_rank)
  data$alpha = data$alpha %>%
    .order_columns(.get_single_chain_cols_ordered(), keep_extra = TRUE, verbose = verbose) %>%
    arrange(rank)
  data$beta = data$beta %>%
    .order_columns(.get_single_chain_cols_ordered(), keep_extra = TRUE, verbose = verbose) %>%
    arrange(rank)

  return(data)
}


.get_paired_cols_ordered = function() {
  c("method", "va", "ja", "cdr3a", "cdr3b", "vb", "jb",
    "alpha_rank", "beta_rank",
    "alpha_readFraction", "beta_readFraction",
    "alpha_readCount", "beta_readCount",
    "is_functional",
    "alpha_nuc", "beta_nuc", "alpha_beta",
    "wi", "wj", "wij", "wa", "wb", "score",
    "r", "ts", "pval", "pval_adj", "loss_a_frac", "loss_b_frac",
    "alpha_readCount_max", "alpha_readCount_median",
    "alpha_sem", "alpha_max_wells",
    "beta_readCount_max", "beta_readCount_median",
    "beta_sem",  "beta_max_wells",
    "alpha_has_stop_codon", "alpha_has_frameshift",
    "beta_has_stop_codon", "beta_has_frameshift",
    "alpha_is_functional", "beta_is_functional"
  )
}

.get_paired_alt_cols_ordered = function() {
  cols = c("paired_status", .get_paired_cols_ordered(), "is_paired_madhype", "is_paired_tshell")
  cols = cols[!cols == "method"]
  return(cols)
}

.get_single_chain_cols_ordered = function() {
  c("id", "rank","aaSeqCDR3", "v", "j",
    "readFraction", "readCount", "n_wells",
    "paired_status", "n_paired",
    "targetSequences",
    "readCount_median", "readCount_max",
    "sem", "max_wells",
    "n_paired_madhype", "n_paired_tshell",
    "has_stop_codon", "has_frameshift", "is_functional",
    "is_paired", "is_paired_tshell",
    "is_paired_madhype")
}

.get_receptor_component = function(dt, method = c("beta_id","components")) {
  method = method[1]
  checkmate::assert_choice(method, choices = c("beta_id","components"))
  if(method == "components") {
    requireNamespace("igraph", quietly = TRUE)
    dt2 <- dt[,c("alpha_id", "beta_id")]
    # Build graph and get connected components
    g     <- igraph::graph_from_data_frame(dt2, directed = FALSE)
    comps <- igraph::components(g)
    # Map component back to original rows
    dt$receptor_rank = comps$membership[dt$alpha_id]
    dt$receptor_id = paste("receptor", dt$receptor_rank, sep = "_")
    #dt$receptor_hash = sapply(dt$receptor_id, function(x) digest::digest(x, algo = "xxhash32"))
  } else if(method == "beta_id") {
    dt$receptor_rank = dt$beta_rank
    dt$receptor_id = paste("receptor", dt$receptor_rank, sep = "_")
    #dt$receptor_hash = sapply(dt$receptor_id, function(x) digest::digest(x, algo = "xxhash32"))
  }
  return(dt)
}

## input is a list with "paired", "alpha", and "beta" slots
## output is a paired data frame
.add_single_chain_data_simple = function(data, verbose=TRUE) {
  sc_cols = c("targetSequences", "readCount", "readCount_max", "readCount_median", "sem", "readFraction", "max_wells", "rank", "id")
  dt_pair = data$paired
  dtA = data$alpha %>% .order_columns(cols = sc_cols, keep_extra = FALSE, verbose = verbose)
  dtB = data$beta %>% .order_columns(cols = sc_cols, keep_extra = FALSE, verbose = verbose)
  colnames(dtA) = paste("alpha_", colnames(dtA), sep = "")
  colnames(dtB) = paste("beta_", colnames(dtB), sep = "")
  colnames(dtA)[1] = "alpha_nuc"
  colnames(dtB)[1] = "beta_nuc"
  if(!is.data.table(dtA)) dtA = as.data.table(dtA)
  if(!is.data.table(dtB)) dtB = as.data.table(dtB)
  if(!is.data.table(dt_pair)) dt_pair = as.data.table(dt_pair)
  dt_join1 = data.table::merge.data.table(x=dt_pair, y=dtA, by = "alpha_nuc", all.x=TRUE, all.y=FALSE)
  dt_join2 = data.table::merge.data.table(x=dt_join1, y=dtB, by = "beta_nuc", all.x=TRUE, all.y=FALSE)

  return(dt_join2)
}

.annotate_paired_single = function(data, by_method=TRUE) {
  ll = .identify_paired_single(data, by_method=by_method)
  data$alpha = ll$alpha
  data$beta = ll$beta
  data$paired_alt = ll$paired_alt ## paired data frame w/ no duplicates
  return(data)
}

# for input of a single sample/experiment, returns alpha, beta, and
# paired (with no duplicates) data frames with methods (madhype and tshell) annotated
.identify_paired_single = function(data, by_method=TRUE) {
  if(.is.DataFrame(data$alpha)) data$alpha = as.data.table(data$alpha)
  if(.is.DataFrame(data$beta)) data$beta = as.data.table(data$beta)
  if(.is.DataFrame(data$paired)) data$paired = as.data.table(data$paired)

  paired_tmp = data$paired %>% remove_duplicates() %>% as.data.table()

  counts_alpha = paired_tmp[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired = N, targetSequences = alpha_nuc)
  counts_beta = paired_tmp[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired = N, targetSequences = beta_nuc)

  if(!by_method) {
    tmp_alpha = data$alpha %>%
      left_join(counts_alpha, by = "targetSequences") %>%
      mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
      mutate(is_paired = n_paired != 0) %>%
      mutate(paired_status = ifelse(is_paired, "paired", "un-paired"))
    tmp_beta = data$beta %>%
      left_join(counts_beta, by = "targetSequences") %>%
      mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
      mutate(is_paired = n_paired != 0) %>%
      mutate(paired_status = ifelse(is_paired, "paired", "un-paired"))
    paired_tmp$paired_status = "paired"
    return(list(alpha = tmp_alpha, beta = tmp_beta, paired_alt = paired_tmp))
  }

  paired_tmp_madhype = data$paired %>% filter(method == "madhype")
  paired_tmp_tshell = data$paired %>% filter(method == "tshell")

  counts_alpha_madhype = paired_tmp_madhype[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired_madhype = N, targetSequences = alpha_nuc)
  counts_beta_madhype = paired_tmp_madhype[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired_madhype = N, targetSequences = beta_nuc)


  counts_alpha_tshell = paired_tmp_tshell[, .N, by = alpha_nuc] %>%
    dplyr::rename(n_paired_tshell = N, targetSequences = alpha_nuc)
  counts_beta_tshell = paired_tmp_tshell[, .N, by = beta_nuc] %>%
    dplyr::rename(n_paired_tshell = N, targetSequences = beta_nuc)

  tmp_alpha = data$alpha %>%
    left_join(counts_alpha_madhype, by = "targetSequences") %>%
    left_join(counts_alpha_tshell, by = "targetSequences") %>%
    left_join(counts_alpha, by = "targetSequences") %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0) %>%
    mutate(n_paired_tshell = ifelse(is.na(n_paired_tshell), 0, n_paired_tshell)) %>%
    mutate(is_paired_tshell = n_paired_tshell != 0) %>%
    mutate(n_paired_madhype = ifelse(is.na(n_paired_madhype), 0, n_paired_madhype)) %>%
    mutate(is_paired_madhype = n_paired_madhype != 0) %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    ))

  tmp_beta = data$beta %>%
    left_join(counts_beta_madhype, by = "targetSequences") %>%
    left_join(counts_beta_tshell, by = "targetSequences") %>%
    left_join(counts_beta, by = "targetSequences") %>%
    mutate(n_paired = ifelse(is.na(n_paired), 0, n_paired)) %>%
    mutate(is_paired = n_paired != 0) %>%
    mutate(n_paired_tshell = ifelse(is.na(n_paired_tshell), 0, n_paired_tshell)) %>%
    mutate(is_paired_tshell = n_paired_tshell != 0) %>%
    mutate(n_paired_madhype = ifelse(is.na(n_paired_madhype), 0, n_paired_madhype)) %>%
    mutate(is_paired_madhype = n_paired_madhype != 0)  %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    ))

  paired_tmp = paired_tmp %>%
    mutate(method = NULL) %>%
    mutate(is_paired_madhype = alpha_beta %in% paired_tmp_madhype$alpha_beta) %>%
    mutate(is_paired_tshell = alpha_beta %in% paired_tmp_tshell$alpha_beta) %>%
    mutate(paired_status = case_when(
      is_paired_madhype & is_paired_tshell ~ "both",
      is_paired_madhype & !is_paired_tshell ~ "MAD-HYPE only",
      is_paired_tshell & !is_paired_madhype ~ "T-SHELL only",
      (!is_paired_tshell) & (!is_paired_madhype) ~ "neither",
      .default = "error"
    )) ## paired data frame with no duplicates, labeled by method

  return(list(alpha = tmp_alpha, beta = tmp_beta, paired_alt = paired_tmp))
}


.identify_non_functional_seqs_single = function(df, remove = FALSE) {
  if(isTRUE(is.null(df))) return(df)
  if(.is.DataFrame(df)) df = as.data.table(df)
  #if("is_functional" %in% colnames(df)) return(df) ## this causes errors
  if("cdr3a" %in% colnames(df) && "cdr3b" %in% colnames(df)) {
    ### for paired data
    df_out = df %>%
      mutate(
        alpha_has_stop_codon = grepl("\\*", cdr3a), # stop codon
        alpha_has_frameshift = grepl("_", cdr3a), # frame shift
        beta_has_stop_codon = grepl("\\*", cdr3b),
        beta_has_frameshift = grepl("_", cdr3b),
      ) %>%
      mutate(
        alpha_is_functional = !(alpha_has_stop_codon | alpha_has_frameshift),
        beta_is_functional = !(beta_has_stop_codon | beta_has_frameshift)
      ) %>%
      mutate(is_functional = alpha_is_functional & beta_is_functional)
  } else {
    ### for pseudobulk data
    df_out = df %>%
      mutate(
        has_stop_codon = grepl("\\*", aaSeqCDR3), # stop codon
        has_frameshift = grepl("_", aaSeqCDR3), # frame shift
      ) %>%
      mutate(is_functional = !(has_stop_codon | has_frameshift))
  }
  if(remove) {
    df_out = df_out %>% filter(is_functional)
  }
  if(!"data.table" %in% class(df_out)) {
    #print("converting to data.table")
    df_out = as.data.table(df_out)
  }
  return(df_out)
}
