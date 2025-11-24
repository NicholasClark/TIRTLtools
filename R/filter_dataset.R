#' Subset a TIRTLseqData object
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' The \code{filter_dataset()} function is used to select a subset of samples from
#' a loaded TIRTLseq dataset and create a new dataset object.
#'
#' @details
#' The function accepts 1) a numeric vector of indices, 2) a character vector of sample names, or
#' 3) a character vector of metadata conditions where each entry is of the form \code{"column==value"}.
#'
#' In the third case, \code{c("cell_type==cd4", "timepoint==tp2")} would, for example, select
#' samples whose \code{cell_type} is \code{cd4} and whose \code{timepoint} is \code{tp2} in
#' the sample metadata.
#'
#' @param data a TIRTLseqData object created by \code{\link{load_tirtlseq}()}
#' @param samples the selected samples. Either 1) a numeric vector of indices,
#' 2) a character vector of sample names, or 3) a character vector of metadata
#' conditions where each entry is of the form \code{"column==value"}.
#'
#' @return
#' A dataset object similar to that created by \code{\link{load_tirtlseq}()},
#' but with only the selected samples.
#'
#' @concept data_wrangling
#' @family data_wrangling
#'
#' @export
#' @examples
#' folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal",
#'   package = "TIRTLtools")
#' sjtrc = load_tirtlseq(folder,
#'   meta_columns = c("marker", "timepoint", "version"), sep = "_",
#'   chain = "paired", verbose = FALSE)
#' p2 = filter_dataset(sjtrc, 1:3) ### by indices
#' p3 = filter_dataset(sjtrc, c("cd8_tp1_v2", "cd8_tp2_v2", "cd8_tp3_v2")) ### by sample names
#' p4 = filter_dataset(sjtrc, "marker==cd4") ### by sample metadata condition
#' p5 = filter_dataset(sjtrc, c("marker==cd4", "timepoint==tp2")) ### by multiple sample metadata conditions
#'

filter_dataset = function(data, samples) {
  if(is.character(samples)) {
    if(sum(!samples %in% data$meta$sample_id) > 0) { ### if samples aren't sample names
      if(sum(!grepl("==", samples)) > 0) stop("One or more samples not found in dataset sample names")
      idx_list = lapply(samples, function(x) { ## if samples are conditions, e.g. "cell_type==cd4"
        vec=strsplit(x, "==")[[1]]
        chk1 = length(vec) != 2
        chk2 = !vec[1] %in% colnames(data$meta)
        if(chk1 || chk2) stop(paste("Error for condition:", x))
        idx = which(data$meta[[ vec[1] ]] == vec[2])
        if(length(idx) == 0) stop("Cannot select zero samples")
        return(idx)
      }) ## list of indices satisfying conditions
      samples = Reduce(intersect, idx_list)
    } else { ### if samples are sample names, get indices
      samples = match(samples, data$meta$sample_id)
    }
  }
  if(is.numeric(samples)) { ### if samples are indices
    data$meta = data$meta[samples,]
    data$data = data$data[samples]
  }
  return(data)
}
