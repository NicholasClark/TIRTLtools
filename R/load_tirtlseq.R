#' Load data from TIRTLseq experiments
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' \code{load_tirtlseq()} loads paired-TCR and pseudo-bulk data from TIRTLseq experiments
#' from a given directory. It can also automatically assemble metadata from filenames.
#'
#' @details
#' The function expects ".tsv" (or ".tsv.gz") files. It looks for files ending in
#' "_pseudobulk_TRA.tsv" (alpha-chain pseudo-bulk), "_pseudobulk_TRB.tsv" (beta-chain pseudo-bulk),
#' and "_TIRTLoutput.tsv" (paired alpha and beta chains).
#'
#' The files in the directory are found using [get_samples()]. If a sample is missing
#' one of its three files, the corresponding data frame is an empty tibble (zero rows)
#' with the expected column names.
#'
#' By default, the function will construct a metadata table with a row for each sample, based
#' on unique strings at the beginning of filenames (before "_TIRTLoutput.tsv" or similar).
#' If the filename contains sample metadata, then it can add multiple columns to the metadata
#' table with this information. For example, if a typical file looks like "cd8_timepoint2_donor1_TIRTLoutput.tsv"
#' and the user supplies \code{c("cell_type", "timepoint", "donor")} for \code{meta_columns} and \code{"_"} for \code{sep},
#' then the metadata table will look like something like this:
#'\preformatted{
#'    sample_id             cell_type   timepoint       donor     label
#'      <chr>               <chr>         <chr>         <chr>     <chr>
#' 1 cd8_timepoint2_donor1    cd8       timepoint2      donor1    cell_type: cd8 | timepoint: timepoint2 | donor: donor1
#' 2 ...
#' 3 cd4_timepoint1_donor3    cd4       timepoint1      donor3    cell_type: cd4 | timepoint: timepoint1 | donor: donor3
#' }
#'
#' @param directory the directory (or a vector of directories) to look in for ".tsv" or ".tsv.gz" files of TIRTLseq data
#' @param chain which TCR chain data to load data for -- "all" chains (alpha, beta, and paired) by default.
#' @param sep (optional) separator in the filename for metadata information ("_" by default)
#' @param meta_columns (optional) a vector of names identifying the metadata contained in the filenames,
#' for example \code{c("marker", "timepoint", "donor")} for files named something like "cd8_timepoint2_donor1 ... .tsv".
#' @param samples (optional) specific sample ids (the part of the filename before "_pseudobulk" or "_TIRTLoutput") to load. Default is NULL (loads all samples in the directory).
#' @param pattern (optional) a pattern used by `grep` to select samples. Default is NULL (loads all samples in the directory).
#' @param clean (optional) a TRUE/FALSE value, whether or not to "clean" the paired data by removing
#' excess pairs for individual alpha and beta chains (default is FALSE).
#' @param remove_nonfunctional whether to remove non-functional TCR chains (default is FALSE)
#' @param process whether to process the raw data, adding columns for alpha and beta readcounts to the
#' paired data frame, adding columns to the pseudobulk data to identify which chains are paired, etc.
#' @param pseudobulk_columns (optional) the columns of the pseudobulk .tsv(.gz) to read. Either
#' a list of columns or one of "auto", "all", or "minimal". "auto" (default) loads all columns except for
#' some redundant ones. "all" loads all columns. "minimal" loads a small number of the most important columns.
#' @param paired_columns (optional) the columns of the paired .tsv(.gz) to read. Either
#' a list of columns or one of "auto", "all", or "minimal". "auto" (default) loads all columns except for
#' some redundant ones. "all" loads all columns. "minimal" loads a small number of the most important columns.
#' @param n_threads (optional) number of CPU threads to use for reading .tsv(.gz) files
#' @param verbose (optional) whether to print the name of each file loaded (default is TRUE).
#' @param stringsAsFactors (optional) read character strings in as factors
#' @param n_max (optional) the maximum number of files to read in -- used mostly for testing purposes (default is Inf, i.e. read all files).
#'
#' @return
#' The function returns a list with the following structure, containing the data from
#' each sample and the sample metadata.
#'
#' \preformatted{
#' your_data_object (list)
#' ├───meta (metadata dataframe)
#' └───data (list)
#'     └───sample_1 (list)
#'         ├───alpha (alpha pseudobulk dataframe)
#'         ├───beta (beta pseudobulk dataframe)
#'         └───paired (paired pseudobulk dataframe)
#'     ...
#'     └───sample_n (list)
#'         ├───alpha (alpha pseudobulk dataframe)
#'         ├───beta (beta pseudobulk dataframe)
#'         └───paired (paired pseudobulk dataframe)
#' }
#' @family data_import
#'
#' @export
#' @examples
#' load_example_data(dataset = "SJTRC_minimal")
#' print(SJTRC_minimal)
#' summary(SJTRC_minimal)
#' SJTRC_minimal$meta ## metadata table
#' names(SJTRC_minimal$data) ## samples
#' SJTRC_minimal$data$cd8_tp1_v2$alpha %>% dplyr::as_tibble() ## alpha chain pseudobulk
#' SJTRC_minimal$data$cd8_tp1_v2$beta %>% dplyr::as_tibble() ## beta chain pseudobulk
#' SJTRC_minimal$data$cd8_tp1_v2$paired %>% dplyr::as_tibble() ## all pairs (rows duplicated if called by both pairing algorithms)
#' SJTRC_minimal$data$cd8_tp1_v2$paired_alt %>% dplyr::as_tibble() ## all pairs (each row unique)
#' SJTRC_minimal$info ## information
#'
load_tirtlseq = function(
    directory,
    chain = c("all","paired","alpha", "beta"),
    sep = "_",
    meta_columns = NULL,
    samples = NULL,
    pattern = NULL,
    clean = FALSE,
    remove_nonfunctional = FALSE,
    process = TRUE,
    pseudobulk_columns = "auto",
    paired_columns = "auto",
    n_threads = data.table::getDTthreads(),
    verbose = TRUE,
    stringsAsFactors = FALSE,
    n_max = Inf) {
  args <- as.list(environment())
  start = Sys.time()
  chain = chain[1]
  compress_strings = FALSE

  checkmate::assert_character(directory, unique = TRUE, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_directory_exists(directory)
  checkmate::assert_character(samples, unique = TRUE, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_string(pattern, null.ok = TRUE, na.ok = FALSE)
  if(!is.null(pattern) && !is.null(samples)) stop("Only one of 'pattern' and 'samples' may be non-null")
  if(!chain %in% c("all","alpha", "beta", "paired")) stop("'chain' must be 'all', 'alpha', 'beta', or 'paired'")
  if("label" %in% meta_columns) stop("'meta_columns' cannot contain a column called 'label'")

  ## the three data frames loaded for each sample: where get_samples() puts the
  ## file path, how to describe it in messages, and which columns to read
  chains = c("alpha", "beta", "paired")
  file_col = c(alpha = "tra_file", beta = "trb_file", paired = "paired_file")
  chain_desc = c(alpha = "Pseudobulk (alpha chain)", beta = "Pseudobulk (beta chain) ", paired = "TIRTLseq (paired chain) ")
  chain_columns = list(alpha = pseudobulk_columns, beta = pseudobulk_columns, paired = paired_columns)

  ### find samples -- one row per sample, with the path to each file (NA if missing)
  if(verbose) message("Loading files from: ", paste(directory, collapse = ", "), "...")
  sample_df = get_samples(directory, pattern = pattern)
  if(verbose) {
    message("Found ", sum(sample_df$has_tra), " alpha chain pseudo-bulk files.")
    message("Found ", sum(sample_df$has_trb), " beta chain pseudo-bulk files.")
    message("Found ", sum(sample_df$has_paired), " paired chain files.")
  }

  ### subset to the samples requested
  if(!is.null(samples)) {
    not_found = setdiff(samples, sample_df$sample)
    if(length(not_found) > 0) {
      warning("No files found for sample(s): ", paste(not_found, collapse = ", "), call. = FALSE)
    }
    sample_df = sample_df[match(setdiff(samples, not_found), sample_df$sample), ]
  }
  if(chain != "all") sample_df = sample_df[!is.na(sample_df[[ file_col[[chain]] ]]), ]
  if(n_max < nrow(sample_df)) {
    message("Loading first ", n_max, " samples:")
    sample_df = sample_df[seq_len(n_max), ]
  }
  sample_ids = sample_df$sample
  if(verbose) message("Loading files for ", length(sample_ids), " samples...")

  ### load the files for each sample
  load_chains = if(chain == "all") chains else chain
  data_final = lapply(seq_along(sample_ids), function(i) {
    if(verbose) message("-- Loading files for sample: ", sample_ids[i], "...")
    obj = lapply(chains, function(chain_tmp) {
      file = sample_df[[ file_col[[chain_tmp]] ]][i]
      columns = chain_columns[[chain_tmp]]
      if(!chain_tmp %in% load_chains) return(.empty_chain_tibble(chain_tmp, columns))
      if(is.na(file)) {
        message("------ File not found -- ", chain_desc[[chain_tmp]])
        return(.empty_chain_tibble(chain_tmp, columns))
      }
      if(verbose) message("---- Loading file -- ", chain_desc[[chain_tmp]], " -- ", basename(file), "...")
      .read_file(file, chain = chain_tmp, cols = columns, n_threads = n_threads,
                 compress_strings = compress_strings, stringsAsFactors = stringsAsFactors)
    }) %>% setNames(chains)
    ## processing needs both pseudobulk data frames
    if( !process || nrow(obj$alpha) == 0 || nrow(obj$beta) == 0 ) return(obj)
    obj = obj %>% .process_TIRTLseq_single(clean = clean, remove_nonfunctional = remove_nonfunctional)
    class(obj) = c("TIRTLseqData", "pairedTCRData")
    return(obj)
  }) %>% setNames(sample_ids)

  n_files = sum(!is.na(unlist(sample_df[ file_col[load_chains] ])))
  if(verbose) message("Loaded ", n_files, " files from ", length(sample_ids), " samples.")

  meta_final = .make_sample_meta(sample_ids, meta_columns = meta_columns, sep = sep)

  duration = Sys.time() - start
  if(verbose) message(paste(round(duration[[1]], 2), units(duration)))
  info = list(data_version = "v1",
              processing_version = "v1",
              package_version = packageVersion("TIRTLtools"),
              call = args)
  out = list(data = data_final, meta = meta_final, info = info)
  class(out) = "TIRTLseqDataSet"
  return(out)
}

## Build the metadata table (one row per sample) from the sample ids, optionally
## splitting each id on `sep` into the columns named in `meta_columns`
.make_sample_meta = function(sample_ids, meta_columns = NULL, sep = "_") {
  if(length(meta_columns) == 0) return(tibble(sample_id = sample_ids, label = sample_ids))
  meta = lapply(sample_ids, function(id) {
    spl = strsplit(id, split = sep)[[1]]
    n = min(length(spl), length(meta_columns))
    c(list(sample_id = id), setNames(as.list(spl[seq_len(n)]), meta_columns[seq_len(n)]))
  }) %>% bind_rows()
  if(nrow(meta) == 0) return(tibble(sample_id = character(), label = character()))
  meta$label = apply(meta[-1], 1, function(row) {
    paste(colnames(meta[-1]), row, sep = ": ", collapse = " | ")
  })
  return(meta)
}

## A zero-row tibble with the columns that would have been read from a file,
## used as a placeholder when a sample has no file for a chain
.empty_chain_tibble = function(chain, cols) {
  all_cols = .get_read_cols(chain, "all")
  read_cols = .get_read_cols(chain, cols)
  if(is.null(read_cols)) read_cols = all_cols
  if(is.null(names(read_cols))) { ## column names only -- look up their types
    types = unname(all_cols[read_cols])
    types[is.na(types)] = "character"
    read_cols = setNames(types, read_cols)
  }
  lapply(read_cols, function(type) {
    tryCatch(vector(type, 0), error = function(e) character(0))
  }) %>% tibble::as_tibble()
}

## The columns (and their types) to read from a pseudobulk ("alpha"/"beta") or
## "paired" file. `cols` is a vector of columns or one of "all", "auto", "minimal".
## Returns NULL (read everything) if `cols` is not recognized.
.get_read_cols = function(chain, cols) {
  if(length(cols) > 1) return(cols)
  if(chain == "paired") {
    all_cols = c(wi = "integer", wj = "integer", wij = "integer", alpha_nuc = "character",
                 beta_nuc = "character", wa = "integer", wb = "integer",
                 alpha_nuc_seq = "character", beta_nuc_seq = "character", alpha_beta = "character",
                 method = "character", r = "double", ts = "double", pval = "double",
                 pval_adj = "double", loss_a_frac = "double", loss_b_frac = "double",
                 score = "double", cdr3a = "character", va = "character", ja = "character",
                 cdr3b = "character", vb = "character", jb = "character")
    auto_drop = c("alpha_nuc_seq", "beta_nuc_seq")
    minimal = c("alpha_nuc", "beta_nuc", "method", "cdr3a", "va", "ja", "cdr3b", "vb", "jb")
  } else {
    all_cols = c(targetSequences = "character", readCount = "integer", v = "character",
                 j = "character", aaSeqCDR3 = "character", n_wells = "integer",
                 readCount_max = "integer", readCount_median = "double", sem = "double",
                 avg = "double",
                 readFraction = "double", max_wells = "integer")
    auto_drop = "avg"
    minimal = c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "sem", "readFraction")
  }
  if(identical(cols, "all")) return(all_cols)
  if(identical(cols, "auto")) return(all_cols[!names(all_cols) %in% auto_drop])
  if(identical(cols, "minimal")) return(all_cols[minimal])
  warning("'cols' must be a list of permissible columns or 'all', 'auto' or 'minimal'. Loading all columns.")
  return(NULL)
}

.read_file = function(file, chain, cols, n_threads, compress_strings, stringsAsFactors) {
  if(!chain %in% c("alpha", "beta", "paired")) {
    message("'chain' must be one of 'alpha', 'beta', or 'paired'. Reading full file.")
    return(data.table::fread(file, nThread = n_threads, stringsAsFactors = stringsAsFactors))
  }
  if(chain %in% c("alpha", "beta")) return(.read_pseudobulk(file, cols, n_threads, compress_strings, stringsAsFactors))
  if(chain %in% c("paired")) return(.read_paired(file, cols, n_threads, compress_strings, stringsAsFactors))
}

.read_paired = function(file, cols, n_threads, compress_strings, stringsAsFactors) {
  read_cols = .get_read_cols("paired", cols)
  df = data.table::fread(file, nThread = n_threads, select = read_cols, stringsAsFactors = stringsAsFactors)
  if(compress_strings) {
    df = DataFrame(df)
    df$alpha_nuc = DNAStringSet(df$alpha_nuc)
    df$beta_nuc = DNAStringSet(df$beta_nuc)
    df$cdr3a = BStringSet(df$cdr3a)
    df$cdr3b = BStringSet(df$cdr3b)
    if("alpha_nuc_seq" %in% colnames(df)) df$alpha_nuc_seq = DNAStringSet(df$alpha_nuc_seq)
    if("beta_nuc_seq" %in% colnames(df)) df$beta_nuc_seq = DNAStringSet(df$beta_nuc_seq)
  }
  return(df)
}

.read_pseudobulk = function(file, cols, n_threads, compress_strings, stringsAsFactors) {
  read_cols = .get_read_cols("pseudobulk", cols)
  df = data.table::fread(file, nThread = n_threads, select = read_cols, stringsAsFactors = stringsAsFactors)
  if(compress_strings) {
    df = DataFrame(df)
    df$targetSequences = DNAStringSet(df$targetSequences)
    df$aaSeqCDR3 = BStringSet(df$aaSeqCDR3)
  }
  return(df)
}
