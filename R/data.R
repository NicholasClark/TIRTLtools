#' Substitution penalty matrix for TCRdist amino acids and V-segments
#'
#' @format A symmetric matrix with substitution penalties for pairs of amino acids or
#' V-segments for TCRdist calculation. The amino acids/V-segments corresponding to each
#' row/column are listed in the "params" table.
#'
#' @concept data
"submat"

#' Table of permissible amino acids and V-segments for TCRdist
#'
#' @format A data frame with two columns. The first column is called "feature" and
#' contains TCR V-segments (e.g. "TRAV12-1*02") and one-letter abbreviations for amino acids.
#' The second column is called "value" and contains the row/column index in the
#' substitution matrix corresponding to that amino acid or V-segment. Note that the
#' first row/column is 0 (instead of 1) because the TCRdist code is in Python, which
#' is 0-indexed.
#' @concept data
"params"
