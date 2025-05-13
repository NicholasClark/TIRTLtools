

## c("dataset", "patient", "timepoint", "cell_type")
load_tirtlseq = function(directory, chain = c("paired","alpha", "beta"), sep = "_", meta_columns = NULL, verbose = TRUE, n_max = Inf) {
  chain = chain[1]
  if(!chain %in% c("alpha", "beta", "paired")) stop("'chain' must be 'alpha' or 'beta', or 'paired'")
  if("label" %in% meta_columns) stop("'meta_columns' cannot contain a column called 'label'")
  if(chain == "alpha") ptn = "*TRA.tsv.gz"
  if(chain == "beta") ptn = "*TRB.tsv.gz"
  if(chain == "paired") ptn = "*TIRTLoutput.tsv.gz"
  files = dir(directory, pattern = ptn)
  if(n_max < length(files)) {
    files = files[1:n_max]
    msg = paste("Loading first ", n_max, " files:", "\n", sep = "")
    cat(msg)
  }
  list_tmp = lapply(files, function(ff) {
    if(verbose) {
      msg = paste("-- Loading file -- ", ff, "\n", sep = "")
      cat(msg)
    }
    data.table::fread(file.path(directory, ff))
    })
  meta_tmp = lapply(files, function(ff) {
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
  return(list(data = list_tmp, meta = meta_tmp))
}
