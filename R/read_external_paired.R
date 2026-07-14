#' Read and process single-cell paired-chain TCR-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' 
#' This function reads and processes paired TCR-sequencing data from a non-TIRTLseq assay.
#' Currently 10X and Parse Biosciences data are supported.
#'
#' @details
#' Supported data types:
#' - \code{"10X"} - "filtered_contig_annotations.csv" or "all_contig_annotations.csv" outputs from Cell Ranger
#' (https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/outputs/cr-5p-outputs-overview-vdj)
#' - \code{"ParseBio"} - "tcr_annotation_airr.tsv" output from Trailmaker
#'
#' @param path the path to the data file
#' @param format the format of the data, either \code{"10X"} or \code{"ParseBio"} or \code{"auto"}. If \code{"auto"} (default), the function will try to decipher which technology the data file was created by.
#' @param id_cols a vector of column names to be used to define a clone (e.g. \code{c("va","vb","ja","jb","alpha_nuc", "beta_nuc", "cdr3a", "cdr3b")}). This can be produced automatically using the \code{\link{make_tcr_schema}()} function.
#' @param multi If \code{FALSE} (default), select only the best two alpha chains for each beta chain
#' when processing the data and creating a data frame with paired TCRs. If \code{TRUE}, keep all alphas.
#' @param separate_rows If \code{TRUE}, when there are multiple alpha chains
#' paired with one beta chain, put each pair in a separate row in the output data frame. If
#' \code{FALSE}, add second alpha chain in extra columns on the same row.
#' @param productive_only If \code{TRUE}, keep only "productive" chains
#' @param verbose If \code{TRUE}, print messages.
#'
#' @return A list containing the following slots:
#' - df_pairs_complete - (data frame) paired receptors - one row for each receptor
#' (excluding those missing an alpha or beta chain)
#' - df_pairs - (data frame) with paired receptors - one row for each receptor
#' (including those missing an alpha or beta chain)
#' - df_pairs_long - (data frame) with paired receptors - one row for each cell
#' (including those missing an alpha or beta chain)
#' - df_pairs_long_complete - (data frame) with paired receptors - one row for each cell
#' (excluding those missing an alpha or beta chain)
#' - df_raw - (data frame) un-edited input data
#' - chain_df - (data frame) summary of total number of each chain in input data
#' - barcode_df - (data frame) summary of number of chains found in each cell
#' - id_cols - (character vector) columns used as IDs to uniquely define receptor pairs
#' - n_cells_total -(integer) total number of cells
#' - n_cells_complete - (integer) number of cells with both chains
#'
#' @family data_processing

read_external_paired = function(
    path,
    format = c("auto","10X", "ParseBio"),
    id_cols = make_tcr_schema(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha = FALSE),
    multi = FALSE,
    separate_rows = TRUE,
    productive_only = TRUE,
    verbose = TRUE) {
  format = format[1]
  df_orig = fread(path)
  alpha_cols = get_tcr_schema_single_chain(id_cols, chain ="alpha")
  beta_cols = get_tcr_schema_single_chain(id_cols, chain ="beta")
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
  df_summ_barcodes = df %>%
    summarize(n_chains = n(),
              n_alphas = sum(chain == "TRA"),
              n_betas = sum(chain == "TRB"),
              n_reads = sum(reads),
              n_umis = sum(umis),
              .by = barcode) %>%
    mutate(complete = (n_alphas > 0) & (n_betas > 0)) %>%
    arrange(desc(n_chains))
  n_cells_complete = sum(df_summ_barcodes$complete)

  complete_barcodes = df_summ_barcodes %>% filter(complete) %>% extract2("barcode")
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
    dplyr::rename(cdr3a = cdr3_aa, alpha_nuc = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis) %>%
    mutate(alpha_receptor_id = paste(!!!syms(alpha_cols), sep = "|")) %>%
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
    dplyr::rename(cdr3b = cdr3_aa, beta_nuc = cdr3_nt, vb = v_gene, jb = j_gene, reads_beta = reads, umis_beta = umis) %>%
    mutate(beta_receptor_id = paste(!!!syms(beta_cols), sep = "|")) %>%
    ungroup()
  df_beta = df_beta_all %>% filter(beta_id == 1)

  # df_beta_multi = df_beta_all %>% filter(beta_id > 1)
  # df_beta_test = df_beta_all %>% filter(n_betas > 1) %>% arrange(barcode)
  # df_alpha_test = df_alpha_all %>% filter(n_alphas > 1) %>% arrange(barcode)

  ## join dfs to make df with top beta and top two alphas for each cell barcode

  ## make df with second alpha for each cell barcode
  if(!separate_rows) {
    df_alpha2 = df_alpha_all %>%
      filter(alpha_id == 2) %>% ## use only second alpha chain
      select(-alpha_id) %>%
      ungroup() %>%
      dplyr::rename(cdr3a2 = cdr3_aa, alpha2_nuc = cdr3_nt, va2 = v_gene, ja2 = j_gene, reads_alpha2 = reads, umis_alpha2 = umis)
    df_pairs_long = df_beta %>%
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
        dplyr::rename(cdr3a = cdr3_aa, alpha_nuc = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis) %>%
        mutate(alpha_receptor_id = paste(!!!syms(alpha_cols), sep = "|"))
      df_alpha12 = bind_rows(df_alpha1, df_alpha2)
      df_pairs_long = df_beta %>%
        full_join(df_alpha12, by=c("barcode")) %>%
        #inner_join(df_alpha12, by=c("barcode")) %>% ## get rid of un-paired chains
        ungroup() %>%
        add_count(barcode, name = "n_in_barcode")
    } else { ## if multi
      df_join = df_alpha_all %>%
        dplyr::rename(cdr3a = cdr3_aa, alpha_nuc = cdr3_nt, va = v_gene, ja = j_gene, reads_alpha = reads, umis_alpha = umis)
      df_pairs_long = df_beta_all %>%
        full_join(df_join, by=c("barcode"))  %>%
        add_count(barcode, name = "n_in_barcode")
    }
  }

  df_pairs_long = df_pairs_long %>%
    mutate(has_beta = !is.na(beta_nuc), has_alpha = !is.na(alpha_nuc)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    #dplyr::rename(alpha_nuc = cdr3a_nt, beta_nuc = cdr3b_nt) %>%
    mutate(receptor_id = paste(!!!syms(id_cols), sep = "|"),
           alpha_receptor_id = paste(!!!syms(alpha_cols), sep = "|"),
           beta_receptor_id = paste(!!!syms(beta_cols), sep = "|")
           ) %>%
    add_count(receptor_id, name = "n_receptor") %>%
    arrange(desc(n_receptor), receptor_id, desc(umis_beta), desc(reads_beta), desc(umis_alpha), desc(reads_alpha))

  if(!separate_rows) df_pairs_long = df_pairs_long %>% dplyr::rename(alpha2_nuc = cdr3a2_nt)

  df_pairs_long_complete = df_pairs_long %>% filter(has_alpha_and_beta)

  # id_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "alpha_nuc", "beta_nuc")
  # if(!separate_rows) id_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "va2", "ja2", "cdr3a2", "alpha_nuc", "beta_nuc", "alpha2_nuc")

  ## df_pairs_long is all pairs of receptors in all cells
  ## df_pairs is a cell count for unique pairs of receptors as defined by "id_cols"
  df_pairs = df_pairs_long %>%
    #mutate(alpha_beta = paste(cdr3a_nt, cdr3b_nt, sep = "_")) %>%
    #add_count(barcode, name = "n_in_barcode") %>%
    summarize(n_cells = sum(1/n_in_barcode), ## number of cells, but dividing if multiple alphas
              n_in_barcode = max(n_in_barcode), ## number of cells not dividing if multiple alphas
              n_cells_in_barcode = max(n_receptor),
              total_umis_alpha = sum(umis_alpha, na.rm = TRUE), total_reads_alpha = sum(reads_alpha, na.rm = TRUE),
              total_umis_beta = sum(umis_beta, na.rm = TRUE), total_reads_beta = sum(reads_beta, na.rm = TRUE),
              .by = c(!!!syms(id_cols),receptor_id, alpha_receptor_id, beta_receptor_id)) %>%
    #add_count(alpha_beta, name = "n_cells") %>%
    mutate(has_alpha = !is.na(alpha_nuc), has_beta = !is.na(beta_nuc)) %>%
    mutate(has_alpha_and_beta = has_alpha & has_beta) %>%
    ## order rows by number of cells, then by beta umis/reads if tied, then if still tied by alpha umis/reads, then by receptor_id
    arrange(desc(n_cells), desc(total_umis_beta), desc(total_reads_beta), desc(total_umis_alpha), desc(total_reads_alpha), receptor_id)

  if(!separate_rows) {
    df_pairs = df_pairs %>% mutate(has_alpha2 = !is.na(alpha2_nuc))
  }
  df_pairs$method = format
  df_pairs_complete = df_pairs %>% filter(has_alpha_and_beta)
  #df_complete <- as.data.table(df_clean)[order(-n_cells),][!duplicated(paste0(beta_nuc,"_",alpha_nuc)),][!is.na(beta_nuc)&!is.na(alpha_nuc),]

  ## summary of how often each complete TCR is found vs. just alpha or just beta
  df_rec = df_pairs_complete %>%
    select(receptor_id, alpha_receptor_id, beta_receptor_id) %>%
    mutate(
      n_cells_with_both = sapply(receptor_id, function(x) sum(df_pairs_long$receptor_id == x)),
      n_cells_with_alpha = sapply(alpha_receptor_id, function(x) sum(df_pairs_long$alpha_receptor_id == x)),
      n_cells_with_beta = sapply(beta_receptor_id, function(x) sum(df_pairs_long$beta_receptor_id == x))
      )

  out = list(df_pairs_complete = df_pairs_complete,    ## (data frame) paired receptors - one row for each receptor (excluding those missing an alpha or beta chain)
             df_pairs = df_pairs,                   ## (data frame) with paired receptors - one row for each receptor (including those missing an alpha or beta chain)
             df_pairs_long = df_pairs_long,              ## (data frame) with paired receptors - one row for each cell (including those missing an alpha or beta chain)
             df_pairs_long_complete = df_pairs_long_complete,     ## (data frame) with paired receptors - one row for each cell (excluding those missing an alpha or beta chain)
             df_raw = df_orig,                     ## (data frame) un-edited input data
             n_cells_total = n_cells_total,        ## (integer) total number of cells
             n_cells_complete = n_cells_complete,  ## (integer) number of cells with both chains
             chain_df = chain_df,                  ## (data frame) summary of total number of each chain in input data
             barcode_df = df_summ_barcodes,        ## (data frame) summary of number of chains found in each cell
             id_cols = id_cols, ## (character vector) columns used as IDs to uniquely define receptor pairs
             n_cells_total = n_cells_total,        ## (integer) total number of cells
             n_cells_complete = n_cells_complete  ## (integer) number of cells with both chains
             )

  return(out)
}

.infer_format = function(df, verbose = TRUE) {
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


.get_necessary_cols_paired = function(format = c("10X", "ParseBio")) {
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
  } else if(format == "10X") {
    ## note: this is for filtered_contig_annotations.csv or all_contig_annotations.csv
    ## "filtered" contig annotations are all of those that are TRUE for "is_cell", "high confidence", "full_length", and "productive"
    ## other 10X/CellRanger output files not included yet:
    full_cols = c("barcode", "is_cell", "contig_id", "high_confidence", "length",
                  "chain", "v_gene", "d_gene", "j_gene", "c_gene", "full_length",
                  "productive", "fwr1", "fwr1_nt", "cdr1", "cdr1_nt", "fwr2", "fwr2_nt",
                  "cdr2", "cdr2_nt", "fwr3", "fwr3_nt", "cdr3", "cdr3_nt", "fwr4",
                  "fwr4_nt", "reads", "umis", "raw_clonotype_id", "raw_consensus_id",
                  "exact_subclonotype_id") ## for reference
    necessary_cols = c("barcode", "chain", "productive", "v_gene", "j_gene", "cdr3_nt", "cdr3", "reads", "umis")
  } else {
    stop("'format' is not recognized: must be either '10X' or 'ParseBio'")
  }
  return(necessary_cols)
}


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
