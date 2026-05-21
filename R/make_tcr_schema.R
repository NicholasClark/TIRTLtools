make_tcr_schema = function(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha = FALSE) {
  cols = c()
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
  return(cols)
}
