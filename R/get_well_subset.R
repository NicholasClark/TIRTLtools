#' Get well names from numerical rows and columns
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' @family well
get_well_subset = function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F)) %>% as.vector()
}
