#' Get well names from top left and bottom right wells
#'
#' `r lifecycle::badge('experimental')`
#'
#' @family well
get_wells_from_edges = function(top_left, bottom_right, return_type = c("wells", "rows_and_columns")){
  return_type = return_type[1]
  row_start = substr(top_left, 1, 1)
  col_start = as.integer( substr(top_left, 2, nchar(top_left)) )
  row_end = substr(bottom_right, 1, 1)
  col_end = as.integer( substr(bottom_right, 2, nchar(bottom_right)) )
  check1 = col_start <= col_end
  if(!check1) stop("top_left column has to be before bottom_right column")
  cols = col_start:col_end
  row_start = match(row_start, LETTERS)
  row_end = match(row_end, LETTERS)
  check2 = row_start <= row_end
  if(!check2) stop("top_left row has to be before bottom_right row")
  rows = row_start:row_end
  row_letters = LETTERS[rows]
  if(return_type == "wells") {
    out = as.vector(unlist(sapply(row_letters,function(x)paste(x,cols,sep=""),simplify = F)))
  } else {
    out = list(rows = row_letters, columns = cols)
  }
  return(out)
}
