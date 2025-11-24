#' Get well names from numerical rows and columns
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function returns a vector of well names from the specified rows and columns.
#'
#' @param row_range a vector of row indices, from 1 (A) to 16 (P)
#' @param col_range a vector of column indices, from 1 to 24
#'
#' @returns A vector of well names covering the specified rows and columns (e.g. "A1", "A2", ... "P24")
#'
#' @family well
#'
#' @examples
#' get_well_subset(1:16, 1:12) ## left half of plate
#'
#' get_well_subset(1:16, 13:24) ## right half of plate
#'
get_well_subset = function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F)) %>% as.vector()
}
