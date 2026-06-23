combine_TIRTLseqData_old = function(data1, data2) {
  check_TIRTLseqData(data1)
  check_TIRTLseqData(data2)
  names1 = data1$meta$sample_id
  names2 = data2$meta$sample_id
  check = length(intersect(names1, names2)) == 0 ## no intersection between names
  if(!check) {
    msg = "Sample names (<data>$meta$sample_id) must not overlap between two data objects"
    stop(msg)
  } else {
    new_meta = bind_rows(data1$meta, data2$meta)
    new_names = new_meta$sample_id
    new_obj = list(data = c(data1$data, data2$data), meta = new_meta)
    names(new_obj$data) = new_names
    class(new_obj) = "TIRTLseqDataSet"
  }
  return(new_obj)
}

combine_TIRTLseqData = function(...) {
  args <- list(...)
  for(x in args) check_TIRTLseqData(x)
  meta = bind_rows(lapply(args, function(x) x$meta))
  tab = table(meta$sample_id)
  tab_multiple = tab[tab>1]
  check = length(tab_multiple) == 0 ## no intersection between names
  if(!check) {
    msg = "Sample names (<data>$meta$sample_id) must not overlap between two data objects"
    stop(msg)
  } else {
    new_obj = list(data = do.call(c, lapply(args, function(x) x$data)), meta = meta)
    new_obj$data = new_obj$data[new_obj$meta$sample_id]
    class(new_obj) = "TIRTLseqDataSet"
  }
  return(new_obj)
}
