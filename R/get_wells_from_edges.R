#' Get well names from the top left and bottom right wells
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function returns a vector of well names in a rectangle containing the top
#' left well and the bottom right well.
#'
#' @param top_left the well in the top left of the rectangle (e.g. "A1")
#' @param bottom_right the well in the top left of the rectangle (e.g. "H12")
#'
#' @returns If return_type = "wells" (default), returns a vector of well names
#' covering the specified rows and columns (e.g. "A1", "A2", ... "H12"). If
#' return_type = "rows_and_columns", returns a list with two vectors
#'  - rows - a vector of the rows (e.g. "A", "B", ..., "H")
#'  - columns - a vector of the columns (e.g. 1, 2, 3, ..., 24)
#'
#' @family well
#'
#' @examples
#' get_wells_from_edges("A13", "P24") ## right half of plate
#'
#' get_wells_from_edges("A1", "H12") ## top left half of plate
#'
#' get_wells_from_edges("A13", "P24", return_type = "rows_and_columns")
#'
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
