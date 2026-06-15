#' Run data processing functions on a TIRTLseqDataSet object
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function runs annotation and data cleaning functions on a
#' TIRTLseqDataSet object. Specifically, it calls the functions \code{add_single_chain_data()},
#' \code{identify_paired()}, \code{identify_non_functional_seqs()}, and (optionally) \code{clean_pairs()}.
#'
#' @param data a TIRTLseqDataSet object
#' @param clean (optional) a TRUE/FALSE value, whether or not to "clean" the paired data by removing
#' excess pairs for individual alpha and beta chains (default is FALSE).
#' @param remove_nonfunctional whether to remove non-functional TCR chains (default is FALSE)
#'
#' @return a TIRTLseqDataSet object with annotated and (optionally) cleaned data
#'
#' @family data_processing
#'
process_TIRTLseq = function(data, clean = FALSE, remove_nonfunctional = FALSE) {
  sample_names = names(data$data)
  for(i in 1:length(data$data)) {
    data$data[[i]] = .process_TIRTLseq_single(data$data[[i]], clean = clean, remove_nonfunctional = remove_nonfunctional)
  }
  names(data$data) = sample_names
  return(data)
}

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
