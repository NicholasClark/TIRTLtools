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
