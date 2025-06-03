#' @title
#' Load TIRTLseq data
#'
#' @description
#' \code{load_tirtlseq()} loads TIRTLseq data from a given directory. It can also automatically
#' assemble metadata from the filenames.
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
#' @param directory the directory to look in for ".tsv" or ".tsv.gz" files
#' @param chain which chain data to load -- all chains (alpha, beta, and paired) by default.
#' @param sep (optional) separator in the filename for metadata information ("_" by default)
#' @param meta_columns (optional) a vector of identifying the metadata contained in filenames,
#' for example \code{c("cell_type", "timepoint", "donor")} for files named similar to "cd8_timepoint2_donor1_TIRTLoutput.tsv".
#' @param verbose whether to print the name of each file loaded (default is TRUE).
#' @param n_max the maximum number of files to read in -- used mostly for testing purposes (default is Inf, i.e. no maximum).
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
#' @seealso \code{\link{filter_dataset}()}
#'
#' @export
#' @examples
#' # example code
#' # paired = load_tirtlseq("path_to/your_directory", sep = "_", meta_columns = c("cell_type", "timepoint"))
#'
#'
load_tirtlseq = function(directory, chain = c("all","paired","alpha", "beta"), sep = "_", meta_columns = NULL, verbose = TRUE, n_max = Inf) {
  tictoc::tic()
  chain = chain[1]
  if(!chain %in% c("all","alpha", "beta", "paired")) stop("'chain' must be 'all', 'alpha', 'beta', or 'paired'")
  if("label" %in% meta_columns) stop("'meta_columns' cannot contain a column called 'label'")

  msg = paste("Loading files from: ", directory, "...", "\n", sep = "")
  cat(msg)
  # if(chain == "alpha") ptn = "*TRA.tsv.gz"
  # if(chain == "beta") ptn = "*TRB.tsv.gz"
  # if(chain == "paired") ptn = "*TIRTLoutput.tsv.gz"
  all_files = dir(directory)
  file_exts = all_files %>% tools::file_ext() %>% unique()

  if(length(file_exts) > 1) stop("All files need to be either .tsv or .tsv.gz")
  if(file_exts == "tsv") {
    alpha_post = "_pseudobulk_TRA.tsv"
    beta_post = "_pseudobulk_TRB.tsv"
    paired_post = "_TIRTLoutput.tsv"
  } else {
    alpha_post = "_pseudobulk_TRA.tsv.gz"
    beta_post = "_pseudobulk_TRB.tsv.gz"
    paired_post = "_TIRTLoutput.tsv.gz"
  }


  files_alpha = dir(directory, pattern = alpha_post )
  files_beta = dir(directory, pattern = beta_post )
  files_paired = dir(directory, pattern = paired_post )

  files_pre_paired = gsub(paired_post, "",files_paired, fixed = TRUE)
  files_pre_alpha = gsub(alpha_post, "",files_alpha, fixed = TRUE)
  files_pre_beta = gsub(beta_post, "",files_beta, fixed = TRUE)

  files_pre = case_when(
    chain == "all" ~ unique(c(files_pre_paired, files_pre_alpha, files_pre_beta)),
    chain == "alpha" ~ files_pre_alpha,
    chain == "beta" ~ files_pre_beta,
    chain == "paired" ~ files_pre_paired
  )

  if(n_max < length(files_pre)) {
    files_pre = files_pre[1:n_max]
    msg = paste("Loading first ", n_max, " files:", "\n", sep = "")
    cat(msg)
  }

  file_counter = 0
  ### loop over each experiment to load files
  list_tmp = lapply(files_pre, function(ff_pre) {
    if(verbose) {
      msg = paste("-- Loading files for sample: ", ff_pre, "...",  "\n", sep = "")
      cat(msg)
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
      file_long = file.path(directory, ff)
      if(chain %in% c("alpha", "beta", "paired") && chain != chain_tmp) return(list())
      if(file.exists(file_long)) {
        if(verbose) {
          msg = paste("---- Loading file -- ", desc, " -- ", ff, "...", "\n", sep = "")
          cat(msg)
          df_tmp =  data.table::fread(file_long)# %>% dtplyr::lazy_dt()
          df_tmp = identify_non_functional_seqs(df_tmp)
          file_counter <<- file_counter + 1
        }
        return(df_tmp)
      } else {
        msg = paste("------ File not found -- ", desc, " -- ", ff, "\n", sep = "")
        cat(msg)
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
      meta_tmp$label = df_tmp$sample_id
    }

  }
  msg = paste("Loaded ", file_counter, " files from ", length(list_tmp), " samples.", "\n", sep = "")
  cat(msg)
  tictoc::toc()
  return(list(data = list_tmp, meta = meta_tmp))
}
