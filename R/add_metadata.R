#' Add metadata to a TIRTLseqData object
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @param obj a TIRTLseqData object
#' @param meta_columns (optional) a vector of names for metadata for each field in the sample names.
#' @param sep (optional) the character separating fields in the sample names of the data
#' For example \code{c("marker", "timepoint", "donor")} for samples named something
#' like "cd8_timepoint2_donor1".
#'
#'
#' @family data_wrangling
#'
add_metadata = function(obj, meta_columns = NULL, sep="_") {
  samples = names(obj$data)
  meta_tmp = lapply(1:length(samples), function(i) {
    ff = samples[i]
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


  meta_tmp$label = apply(meta_tmp[-1], 1, function(row) {
    paste(colnames(meta_tmp[-1]), row, sep = ": ", collapse = " | ")
  })
  obj$meta = meta_tmp
  return(obj)
}
