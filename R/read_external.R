## example:
# path = "data/BM03_processed-data/neb/N-F5L/N-F5L_TRA.tsv"
# path = "data/BM03_processed-data/cellecta_RNA/C2-R-F5L/C2-R-F5L_TRA.tsv"
# path = "data/BM03_processed-data/takara/T-F5L/T-F5L_TRA.tsv"
# path = "data/BM03_processed-data/qiaseq/Q-combined-F5L_combined_run2/Q-combined-F5L.clones_TRA.tsv"
read_external_bulk = function(path) {
  df_orig = fread(path)
  cols = .get_necessary_cols_bulk()
  .check_cols(df_orig, cols) ## errors if it doesn't have the needed columns
  df = df_orig %>%
    mutate(v = gsub("\\*.*", "", allVHitsWithScore),
           j = gsub("\\*.*", "", allJHitsWithScore)) %>%
    select(targetSequences, readCount, v, j, aaSeqCDR3, readFraction, everything())
  return(list(df = df, df_raw = df_orig))
}

## multi -- if FALSE, only take one beta and best two alpha, else take all pairs
## separate_rows -- if TRUE, report each pair on a separate row, instead of extra columns for second alpha
read_external_paired = function(
    path,
    format = c("auto","10x", "ParseBio"),
    id_cols = make_tcr_schema(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha = FALSE),
    multi = FALSE,
    separate_rows = TRUE,
    productive_only = TRUE,
    verbose = TRUE) {
  format = format[1]
  df_orig = fread(path)
  if(format == "auto") format = .infer_format(df_orig)
  cols = .get_necessary_cols_paired(format)
  .check_cols(df_orig, cols) ## errors if it doesn't have the needed columns
  if(format=="ParseBio") df_orig = df_orig %>% mutate(productive = .convert_char_to_boolean(productive), rev_comp = .convert_char_to_boolean(rev_comp))
  df = rename_columns(df_orig, format=format, rename_df = get_names_df_single_chain(), verbose = TRUE)
  if(productive_only) df = df %>% filter(productive)

  n_cells_total = length(unique(df$barcode))
  chain_df = df %>% summarize(n = n(), .by=chain) %>% arrange(chain)
  NA_chains = df %>% filter(is.na(chain)) ## should be nothing with NA for chain
  df = df %>% filter(!is.na(chain))
  df_summ = df %>%
    summarize(n_chains = n(), n_alphas = sum(chain == "TRA"), n_betas = sum(chain == "TRB"), .by = barcode) %>%
    mutate(complete = (n_alphas > 0) & (n_betas > 0)) %>%
    arrange(desc(n_chains))
  n_cells_complete = sum(df_summ$complete)

  complete_barcodes = df_summ %>% filter(complete) %>% extract2("barcode")
  df_complete_tcr_only = df %>% filter(barcode %in% complete_barcodes)

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

  # df_alpha_multi = df_alpha_all %>%
  #   filter(alpha_id > 2) %>% ## second alpha chain
  #   dplyr::rename(cdr3a2 = cdr3_aa, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis) %>%
  #   select(-alpha_id) %>%
  #   ungroup()

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
      dplyr::rename(cdr3a2 = cdr3_aa, cdr3a2_nt = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis)
    df_all1 = df_beta %>%
      full_join(df_alpha1, by=c("barcode")) %>%
      full_join(df_alpha2, by=c("barcode")) %>%
      ungroup() %>%
      add_count(barcode, name = "n_in_barcode")
  } else { ## if separate_rows
    if(!multi) {
      df_alpha2 = df_alpha_all %>%
        filter(alpha_id == 2) %>% ## use only second alpha chain
        select(-alpha_id) %>%
        ungroup() %>%
        dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_alpha12 = bind_rows(df_alpha1, df_alpha2)
      df_all1 = df_beta %>%
        full_join(df_alpha12, by=c("barcode")) %>%
        #inner_join(df_alpha12, by=c("barcode")) %>% ## get rid of un-paired chains
        ungroup() %>%
        add_count(barcode, name = "n_in_barcode")
    } else { ## if multi
      df_join = df_alpha_all %>%
        dplyr::rename(cdr3a = cdr3_aa, cdr3a_nt = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_all1 = df_beta_all %>%
        full_join(df_join, by=c("barcode"))  %>%
        add_count(barcode, name = "n_in_barcode")
    }
  }

  df_all1 = df_all1 %>%
    mutate(has_beta = !is.na(cdr3b_nt), has_alpha = !is.na(cdr3a_nt)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    dplyr::rename(alpha_nuc = cdr3a_nt, beta_nuc = cdr3b_nt) %>%
    mutate(receptor_id = paste(!!!syms(id_cols), sep = "_")) %>%
    add_count(receptor_id, name = "n_receptor") %>%
    arrange(desc(n_receptor), receptor_id, desc(umis_beta), desc(reads_beta), desc(umis_alpha), desc(reads_alpha))

  if(!separate_rows) df_all1 = df_all1 %>% dplyr::rename(alpha2_nuc = cdr3a2_nt)

  df_all2 = df_all1 %>% filter(has_alpha_and_beta)

  # id_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "alpha_nuc", "beta_nuc")
  # if(!separate_rows) id_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "va2", "ja2", "cdr3a2", "alpha_nuc", "beta_nuc", "alpha2_nuc")

  ## df_all1 is all pairs of receptors in all cells
  ## df_summ is a cell count for unique pairs of receptors as defined by "id_cols"
  df_summ = df_all1 %>%
    #mutate(alpha_beta = paste(cdr3a_nt, cdr3b_nt, sep = "_")) %>%
    #add_count(barcode, name = "n_in_barcode") %>%
    summarize(n_cells = sum(1/n_in_barcode),
              total_umis_alpha = sum(umis_alpha, na.rm = TRUE), total_reads_alpha = sum(reads_alpha, na.rm = TRUE),
              total_umis_beta = sum(umis_beta, na.rm = TRUE), total_reads_beta = sum(reads_beta, na.rm = TRUE),
              .by = c(!!!syms(id_cols),receptor_id)) %>%
    #add_count(alpha_beta, name = "n_cells") %>%
    mutate(has_alpha = !is.na(alpha_nuc), has_beta = !is.na(beta_nuc)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    ## order rows by number of cells, then by beta umis/reads if tied, then if still tied by alpha umis/reads, then by receptor_id
    arrange(desc(n_cells), desc(total_umis_beta), desc(total_reads_beta), desc(total_umis_alpha), desc(total_reads_alpha), receptor_id)

  if(!separate_rows) {
    df_summ = df_summ %>% mutate(has_alpha2 = !is.na(alpha2_nuc))
  }
  df_summ_clean = df_summ %>% filter(has_alpha_and_beta)
  #df_complete <- as.data.table(df_clean)[order(-n_cells),][!duplicated(paste0(beta_nuc,"_",alpha_nuc)),][!is.na(beta_nuc)&!is.na(alpha_nuc),]

  out = list(df_pairs_complete = df_summ_clean,    ## df with number of cells for each paired receptor (with both alpha & beta)
             df_pairs = df_summ,                   ## df with number of cells for each paired receptor
             df_pairs_long = df_all1,              ## df with all individual paired receptors
             df_pairs_long_complete = df_all2,     ## df with all individual paired receptors (with both alpha & beta)
             df_raw = df_orig,                     ## un-edited input data
             n_cells_total = n_cells_total,        ## total number of cells
             n_cells_complete = n_cells_complete,  ## number of cells with both chains
             chain_df = chain_df,                  ## summary of total number of each chain in input data
             barcode_df = df_summ,                 ## summary of number of chains found in each cell
             id_cols = id_cols)                    ## columns used as IDs to uniquely define receptor pairs

  return(out)
}

.infer_format = function(df, verbose = TRUE) {
  if("locus" %in% colnames(df)) {
    if(verbose) message("Reading as 'ParseBio'...")
    return("ParseBio")
  }
  if("chain" %in% colnames(df)) {
    if(verbose) message("Reading as '10x'...")
    return("10x")
  }
  if(verbose) message("Format could not be identified, trying to read as '10x'...")
  return("10x")
}

#' list columns needed in each format
.get_necessary_cols_paired = function(format = c("10x", "ParseBio")) {
  format = format[1]
  if(format == "ParseBio") {
    ## note: this is for tcr_annotation_airr.tsv
    ## other ParseBio output files not included yet: barcode_report.tsv, clonotype_frequency.tsv, tcr_contigs.fa
    full_cols = c("sequence_id", "sequence", "locus", "rev_comp", "productive",
                  "v_call", "d_call", "j_call", "c_call", "sequence_alignment",
                  "germline_alignment", "cdr1", "cdr2", "cdr3", "cdr3_aa", "v_cigar",
                  "d_cigar", "j_cigar", "v_identity", "j_identity", "cell_barcode",
                  "complete_vdj", "read_count", "transcript_count") ## for reference
    necessary_cols =c("cell_barcode", "locus", "productive", "v_call", "j_call", "cdr3", "cdr3_aa", "read_count", "transcript_count")
  } else if(format == "10x") {
    ## note: this is for filtered_contig_annotations.csv or all_contig_annotations.csv
    ## "filtered" contig annotations are all of those that are TRUE for "is_cell", "high confidence", "full_length", and "productive"
    ## other 10x/CellRanger output files not included yet:
    full_cols = c("barcode", "is_cell", "contig_id", "high_confidence", "length",
                  "chain", "v_gene", "d_gene", "j_gene", "c_gene", "full_length",
                  "productive", "fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt",
                  "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "cdr3", "cdr3_nt", "fwr4",
                  "fwr4_nt", "reads", "umis", "raw_clonotype_id", "raw_consensus_id",
                  "exact_subclonotype_id") ## for reference
    necessary_cols = c("barcode", "chain", "productive", "v_gene", "j_gene", "cdr3_nt", "cdr3", "reads", "umis")
  } else {
    stop("'format' is not recognized: must be either '10x' or 'ParseBio'")
  }
  return(necessary_cols)
}

#' list columns needed in each format
.get_necessary_cols_bulk = function() {
  full_cols = c("cloneId", "readCount", "readFraction", "uniqueMoleculeCount",
    "uniqueMoleculeFraction", "targetSequences", "targetQualities",
    "allVHitsWithScore", "allDHitsWithScore", "allJHitsWithScore",
    "allCHitsWithScore", "allVAlignments", "allDAlignments", "allJAlignments",
    "allCAlignments", "nSeqCDR3", "minQualCDR3", "aaSeqCDR3", "refPoints",
    "uniqueUMICount")
  necessary_cols = c("targetSequences", "readCount", "readFraction", "uniqueMoleculeCount", "uniqueMoleculeFraction", "aaSeqCDR3")
  ## note: might also want targetQualities and minQualCDR3
  return(necessary_cols)
}

.check_cols = function(df, cols) {
  missing_cols = cols[!cols %in% colnames(df)]
  if(length(missing_cols) > 0) {
    msg = paste("Missing columns:", paste(missing_cols, collapse = ", "))
    stop(msg)
  } else {
    return(invisible(NULL))
  }
}

.convert_char_to_boolean = function(vec) {
  if(class(vec) == "logical") return(vec)
  trues = c("T", "TRUE", "True", "true")
  falses = c("F", "FALSE", "False", "false")
  new_vec = case_when(
    vec %in% trues ~ TRUE,
    vec %in% falses ~ FALSE,
    .default = NA
  )
  return(new_vec)
}
