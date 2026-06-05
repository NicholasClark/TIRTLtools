make_tcr_schema = function(features = c("v", "j", "cdr3_aa", "cdr3_nt"), second_alpha = FALSE, rearrange=TRUE) {
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
  if(rearrange) cols = permitted_cols[permitted_cols %in% cols]
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
