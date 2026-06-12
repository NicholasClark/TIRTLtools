#' Make a "schema" for defining a T-Cell Receptor
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This is a convenience function that simply returns a vector of column names
#' that are used by \code{\link{read_external_paired}()} to define a unique
#' T-Cell receptor.
#'
#' For example, using \code{features = c("v", "cdr3_aa")} will group all chains
#' with the same V-alpha/beta and CDR3-alpha/beta amino acid sequence together, while
#' using \code{features = c("v", "cdr3_aa", "cdr3_nt")} will require that their
#' nucleotide sequence is also matching.
#'
#' @param features some subset of \code{c("v", "j", "cdr3_aa", "cdr3_nt")} that
#' you would like to use to define a unique receptor.
#' @param second_alpha whether to include a second alpha chain (default is FALSE)
#'
#' @return A vector containing
#'
#' @examples
#' make_tcr_schema() ## default: v, j, and cdr3 amino acid and nucleotide sequence
#' make_tcr_schema(c("v", "j")) ## only v and j segments
#'
#' @family data_processing

make_tcr_schema = function(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha = FALSE) {
  cols = c()
  permitted_features = c("v", "j", "cdr3_aa", "cdr3_nt")
  check1 = sum(!features %in% permitted_features) > 0
  if(check1) {
    msg = paste("'features' must be selected from the following:", paste(permitted_features, collapse = ", "))
    stop(msg)
  }
  permitted_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "va2", "ja2", "cdr3a2", "alpha_nuc", "beta_nuc", "alpha2_nuc")
  if("v" %in% features) cols = c(cols, "va", "vb")
  if("j" %in% features) cols = c(cols, "ja", "jb")
  if("cdr3_aa" %in% features) cols = c(cols, "cdr3a", "cdr3b")
  if("cdr3_nt" %in% features) cols = c(cols, "alpha_nuc", "beta_nuc")
  if(second_alpha) {
    if("v" %in% features) cols = c(cols, "va2")
    if("j" %in% features) cols = c(cols, "ja2")
    if("cdr3_aa" %in% features) cols = c(cols, "cdr3a2")
    if("cdr3_nt" %in% features) cols = c(cols, "alpha2_nuc")
  }
  cols = permitted_cols[permitted_cols %in% cols]
  return(cols)
}

get_tcr_schema_single_chain = function(schema, chain) {
  chain = chain[1]
  permitted_cols = c("va", "ja", "cdr3a", "cdr3b", "vb", "jb", "va2", "ja2", "cdr3a2", "alpha_nuc", "beta_nuc", "alpha2_nuc")
  alpha_cols_all = c("va", "ja", "cdr3a", "va2", "ja2", "cdr3a2", "alpha_nuc", "alpha2_nuc")
  beta_cols_all = c("vb", "jb", "cdr3b", "beta_nuc")
  alpha_cols = alpha_cols_all[alpha_cols_all %in% schema]
  beta_cols = beta_cols_all[beta_cols_all %in% schema]
  if(chain == "alpha") return(alpha_cols)
  if(chain == "beta") return(beta_cols)
}
