

## c("dataset", "patient", "timepoint", "cell_type")
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
    alpha_post = "_pseudobulk_TRA\\.tsv"
    beta_post = "_pseudobulk_TRB\\.tsv"
    paired_post = "_TIRTLoutput\\.tsv"
  } else {
    alpha_post = "_pseudobulk_TRA\\.tsv\\.gz"
    beta_post = "_pseudobulk_TRB\\.tsv\\.gz"
    paired_post = "_TIRTLoutput\\.tsv\\.gz"
  }


  files_alpha = dir(directory, pattern = paste("*",alpha_post,sep="") )
  files_beta = dir(directory, pattern = paste("*",beta_post,sep="") )
  files_paired = dir(directory, pattern = paste("*",paired_post,sep="") )

  files_pre_paired = gsub(paired_post, "",files_paired)
  files_pre_alpha = gsub(alpha_post, "",files_alpha)
  files_pre_beta = gsub(beta_post, "",files_beta)

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
          df_tmp =  data.table::fread(file_long)
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
  meta_tmp = lapply(files_pre, function(ff) {
    spl = strsplit(ff, split = sep)[[1]]
    df_tmp = list()
    df_tmp$filename = ff
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
    meta_tmp$label = df_tmp$filename
  }
  msg = paste("Loaded ", file_counter, " files from ", length(list_tmp), " samples.", "\n", sep = "")
  cat(msg)
  tictoc::toc()
  return(list(data = list_tmp, meta = meta_tmp))
}
