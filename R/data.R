#' Substitution penalty matrix for TCRdist amino acids and V-segments
#'
#' @format A symmetric matrix (161 x 161) with substitution penalties for pairs of amino acids or
#' V-segments for TCRdist calculation. The amino acids/V-segments corresponding to each
#' row/column are listed in the "params" table.
#'
#' @concept data
#'
#' @examples
#' dim(submat)
#' submat[1:5,1:5]
#' ## to find the amino acid or v-segment(s) that map to a specific row/column
#' index = 106 ## 106th row/column of matrix
#' params[params$value == index-1,] ## need to subtract 1 because values in params are 0-indexed
#'
"submat"

#' Table of permissible amino acids and V-segments for TCRdist
#'
#' @format A data frame with two columns. The data frame contains 266 features, mapping to 161 row/column indices.
#'
#' The first column is called "feature" and
#' contains TCR V-segments (e.g. "TRAV12-1*02") and one-letter abbreviations for amino acids.
#' The second column is called "value" and contains the row/column index in the
#' substitution matrix corresponding to that amino acid or V-segment. Note that the
#' first row/column is 0 (instead of 1) because the TCRdist code is in Python, which
#' is 0-indexed.
#'
#' @concept data
#'
#' @examples
#' params
#'
"params"
