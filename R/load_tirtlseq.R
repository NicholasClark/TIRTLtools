#' Load data from TIRTLseq experiments
#'
#' @description
#' \code{load_tirtlseq()} loads paired-TCR and pseudo-bulk data from TIRTLseq experiments
#' from a given directory. It can also automatically assemble metadata from filenames.
#'
#' @details
#' The function expects ".tsv" (or ".tsv.gz") files. It looks for files ending in
#' "_pseudobulk_TRA.tsv" (alpha-chain pseudo-bulk), "_pseudobulk_TRB.tsv" (beta-chain pseudo-bulk),
#' and "_TIRTLoutput.tsv" (paired alpha and beta chains).
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
#' @param directory the directory to look in for ".tsv" or ".tsv.gz" files of TIRTLseq data
#' @param chain which TCR chain data to load data for -- "all" chains (alpha, beta, and paired) by default.
#' @param sep (optional) separator in the filename for metadata information ("_" by default)
#' @param meta_columns (optional) a vector of names identifying the metadata contained in the filenames,
#' for example \code{c("marker", "timepoint", "donor")} for files named something like "cd8_timepoint2_donor1 ... .tsv".
#' @param samples (optional) specific sample ids (the part of the filename before "_pseudobulk" or "_TIRTLoutput") to load. Default is NULL (loads all samples in the directory).
#' @param pseudobulk_columns (optional) the columns of the pseudobulk .tsv(.gz) to read. Either
#' a list of columns or one of "auto", "all", or "minimal". "auto" (default) loads all columns except for
#' some redundant ones. "all" loads all columns. "minimal" loads a small number of the most important columns.
#' @param paired_columns (optional) the columns of the paired .tsv(.gz) to read. Either
#' a list of columns or one of "auto", "all", or "minimal". "auto" (default) loads all columns except for
#' some redundant ones. "all" loads all columns. "minimal" loads a small number of the most important columns.
#' @param n_threads (optional) number of CPU threads to use for reading .tsv(.gz) files
#' @param compress_strings (optional) whether to compress nucleotide and amino acid sequences using the Biostrings package.
#' @param verbose (optional) whether to print the name of each file loaded (default is TRUE).
#' @param stringsAsFactors (optional) read character strings in as factors
#' @param n_max (optional) the maximum number of files to read in -- used mostly for testing purposes (default is Inf, i.e. read all files).
#'
#' @return
#' The function returns a list with two objects:
#'
#' \code{$meta} - a metadata table (data frame)
#'
#' \code{$data} - a list with one entry for each sample. Each entry is a list with entries
#' \code{$alpha}, \code{$beta}, and \code{$paired}, which are data frames for the alpha- and beta-chain
#' pseudo-bulk data and the paired data respectively.
#'
#' @family data_wrangling
#'
#' @export
#' @examples
#' # example code
#' # paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
#'
#'
load_tirtlseq = function(
    directory,
    chain = c("all","paired","alpha", "beta"),
    sep = "_",
    meta_columns = NULL,
    samples = NULL,
    pseudobulk_columns = "auto",
    paired_columns = "auto",
    n_threads = data.table::getDTthreads(),
    compress_strings = FALSE,
    verbose = TRUE,
    stringsAsFactors = FALSE,
    n_max = Inf) {
  tictoc::tic()
  chain = chain[1]
  if(!chain %in% c("all","alpha", "beta", "paired")) stop("'chain' must be 'all', 'alpha', 'beta', or 'paired'")
  if("label" %in% meta_columns) stop("'meta_columns' cannot contain a column called 'label'")
  ll = lapply(directory, function(dir_tmp) {
    msg = paste("Loading files from: ", dir_tmp, "...", sep = "")
    message(msg)
    # if(chain == "alpha") ptn = "*TRA.tsv.gz"
    # if(chain == "beta") ptn = "*TRB.tsv.gz"
    # if(chain == "paired") ptn = "*TIRTLoutput.tsv.gz"
    all_files = dir(dir_tmp, pattern = ".*(_pseudobulk_TRA\\.tsv|_pseudobulk_TRB\\.tsv|_TIRTLoutput\\.tsv)")
    file_exts = all_files %>% tools::file_ext() %>% unique()

    #if(length(file_exts) > 1) stop("All files need to be either .tsv or .tsv.gz")
    alpha_post = "_pseudobulk_TRA.tsv"
    beta_post = "_pseudobulk_TRB.tsv"
    paired_post = "_TIRTLoutput.tsv"
    # if(file_exts == "tsv") {
    #   alpha_post = "_pseudobulk_TRA.tsv"
    #   beta_post = "_pseudobulk_TRB.tsv"
    #   paired_post = "_TIRTLoutput.tsv"
    # } else {
    #   alpha_post = "_pseudobulk_TRA.tsv.gz"
    #   beta_post = "_pseudobulk_TRB.tsv.gz"
    #   paired_post = "_TIRTLoutput.tsv.gz"
    # }


    files_alpha = dir(dir_tmp, pattern = alpha_post )
    files_beta = dir(dir_tmp, pattern = beta_post )
    files_paired = dir(dir_tmp, pattern = paired_post )

    files_pre_paired = gsub("_TIRTLoutput\\.tsv.*", "",files_paired)
    files_pre_alpha = gsub("_pseudobulk_TRA\\.tsv.*", "",files_alpha)
    files_pre_beta = gsub("_pseudobulk_TRB\\.tsv.*", "",files_beta)

    files_pre = case_when(
      chain == "all" ~ unique(c(files_pre_paired, files_pre_alpha, files_pre_beta)),
      chain == "alpha" ~ files_pre_alpha,
      chain == "beta" ~ files_pre_beta,
      chain == "paired" ~ files_pre_paired
    )

    if(n_max < length(files_pre)) {
      files_pre = files_pre[1:n_max]
      msg = paste("Loading first ", n_max, " files:", sep = "")
      message(msg)
    }

    if(!is.null(samples)) {
      files_pre = samples
    }

    file_counter = 0
    ### loop over each experiment to load files
    list_tmp = lapply(files_pre, function(ff_pre) {
      if(verbose) {
        msg = paste("-- Loading files for sample: ", ff_pre, "...", sep = "")
        message(msg)
      }
      fa = paste(ff_pre, alpha_post, sep = "") %>% gsub("\\\\", "", .)
      fb = paste(ff_pre, beta_post, sep = "") %>% gsub("\\\\", "", .)
      fp = paste(ff_pre, paired_post, sep = "") %>% gsub("\\\\", "", .)
      obj_names = c("alpha", "beta", "paired")
      obj_desc = c("Pseudobulk (alpha chain)", "Pseudobulk (beta chain) ", "TIRTLseq (paired chain) ")
      ff_all = c(fa, fb, fp)
      ### load files for one experiment
      obj = lapply(1:length(ff_all), function(i) {
        desc = obj_desc[i]
        chain_tmp = obj_names[i]
        ff = ff_all[i]
        if(!file.exists(file.path(dir_tmp, ff))) {
          file_long = paste(file.path(dir_tmp, ff), ".gz", sep = "")
        } else {
          file_long = file.path(dir_tmp, ff)
        }


        if(chain_tmp %in% c("alpha", "beta")) columns = pseudobulk_columns
        if(chain_tmp %in% "paired") columns = paired_columns
        if(chain %in% c("alpha", "beta", "paired") && chain != chain_tmp) return(list())
        if(file.exists(file_long)) {
          if(verbose) {
            msg = paste("---- Loading file -- ", desc, " -- ", ff, "...", sep = "")
            message(msg)
          }
          df_tmp = .read_file(file_long, chain = chain_tmp, cols = columns, n_threads = n_threads, compress_strings = compress_strings, stringsAsFactors = stringsAsFactors)
          #df_tmp =  data.table::fread(file_long, nThread = n_threads) # %>% dtplyr::lazy_dt()
          #if(chain %in% c("alpha", "beta")) setkey(df_tmp, targetSequences)
          #df_tmp = identify_non_functional_seqs(df_tmp)
          file_counter <<- file_counter + 1
          return(df_tmp)
        } else {
          msg = paste("------ File not found -- ", desc, " -- ", ff, sep = "")
          message(msg)
          return(list())
        }
      }) %>% setNames(obj_names)
      return(obj)
      }) %>% setNames(files_pre)
    ### make metadata table
    if(is.null(meta_columns)) {
      meta_tmp = tibble(sample_id = files_pre, label = files_pre)
    } else {
      meta_tmp = lapply(1:length(files_pre), function(i) {
        ff = files_pre[i]
        spl = strsplit(ff, split = sep)[[1]]
        df_tmp = list()
        df_tmp$sample_id = ff
        if(length(meta_columns) > 0) {
          for(i in 1:length(meta_columns)) {
            if(i > length(spl)) next
            col = meta_columns[i]
            df_tmp[[ col ]] = spl[i]
          }
        }
        return(df_tmp)
      }) %>% bind_rows()

      if(length(meta_columns) > 0) {
        meta_tmp$label = apply(meta_tmp[-1], 1, function(row) {
          paste(colnames(meta_tmp[-1]), row, sep = ": ", collapse = " | ")
        })
      } else {
        meta_tmp$label = meta_tmp$sample_id
      }

    }
    msg = paste("Loaded ", file_counter, " files from ", length(list_tmp), " samples.", sep = "")
    message(msg)
    return(list(data = list_tmp, meta = meta_tmp))
  })
  meta_final = lapply(ll, function(x) x$meta) %>% bind_rows()
  data_final = lapply(ll, function(x) x$data) %>% do.call(c, .)
  tictoc::toc()
  return(list(data = data_final, meta = meta_final))
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
  if(length(cols) > 1) {
    read_cols = cols
  } else if(cols == "all") {
    #read_cols = NULL
    read_cols = c(wi = "integer", wj = "integer", wij = "integer", alpha_nuc = "character",
                  beta_nuc = "character", wa = "integer", wb = "integer",
                  alpha_nuc_seq = "character", beta_nuc_seq = "character", alpha_beta = "character",
                  method = "character", r = "double", ts = "double", pval = "double",
                  pval_adj = "double", loss_a_frac = "double", loss_b_frac = "double",
                  score = "double", cdr3a = "character", va = "character", ja = "character",
                  cdr3b = "character", vb = "character", jb = "character")
  } else if(cols == "auto") {
    # read_cols = c("wi", "wj", "wij", "alpha_nuc", "beta_nuc", "wa", "wb",
    #          "alpha_beta", "method", "r", "ts", "pval", "pval_adj",
    #          "loss_a_frac", "loss_b_frac", "score", "cdr3a", "va", "ja", "cdr3b",
    #          "vb", "jb")
    read_cols = c(wi = "integer", wj = "integer", wij = "integer", alpha_nuc = "character",
                  beta_nuc = "character", wa = "integer", wb = "integer", alpha_beta = "character",
                  method = "character", r = "double", ts = "double", pval = "double",
                  pval_adj = "double", loss_a_frac = "double", loss_b_frac = "double",
                  score = "double", cdr3a = "character", va = "character", ja = "character",
                  cdr3b = "character", vb = "character", jb = "character")
  } else if(cols == "minimal") {
    # read_cols = c("alpha_nuc", "beta_nuc", "method",
    #               "cdr3a", "va", "ja", "cdr3b", "vb", "jb")
    read_cols = c(alpha_nuc = "character", beta_nuc = "character",
      method = "character", cdr3a = "character", va = "character", ja = "character",
      cdr3b = "character", vb = "character", jb = "character")
  } else {
    msg = "'cols' must be a list of permissible columns or 'all', 'auto' or 'minimal'. Loading all columns."
    warning(msg)
    read_cols = NULL
  }
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
  if(length(cols) > 1) {
    read_cols = cols
  } else if(cols == "all") {
    #read_cols = NULL
    read_cols = c(targetSequences = "character", readCount = "integer", v = "character",
                  j = "character", aaSeqCDR3 = "character", n_wells = "integer",
                  readCount_max = "integer", readCount_median = "double", sem = "double",
                  avg = "double",
                  readFraction = "double", max_wells = "integer")
  } else if(cols == "auto") {
    # read_cols = c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "n_wells",
    #               "readCount_max", "readCount_median", "sem", "readFraction",
    #               "max_wells")
    read_cols = c(targetSequences = "character", readCount = "integer", v = "character",
                  j = "character", aaSeqCDR3 = "character", n_wells = "integer",
                  readCount_max = "integer", readCount_median = "double", sem = "double",
                  readFraction = "double", max_wells = "integer")
  } else if(cols == "minimal") {
    #read_cols = c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "sem", "readFraction")
    read_cols = c(targetSequences = "character", readCount = "integer", v = "character",
      j = "character", aaSeqCDR3 = "character", sem = "double",
      readFraction = "double")
  } else {
    msg = "'cols' must be a list of permissible columns or 'all', 'auto' or 'minimal'. Loading all columns."
    warning(msg)
    read_cols = NULL
  }
  df = data.table::fread(file, nThread = n_threads, select = read_cols, stringsAsFactors = stringsAsFactors)
  if(compress_strings) {
    df = DataFrame(df)
    df$targetSequences = DNAStringSet(df$targetSequences)
    df$aaSeqCDR3 = BStringSet(df$aaSeqCDR3)
  }
  return(df)
}
