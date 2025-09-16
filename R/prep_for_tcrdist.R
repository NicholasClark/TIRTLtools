### prep data for tcrdist python function
prep_for_tcrdist = function(df, params=NULL) {
  if(is.null(df)) return(df)
  if(is.null(params)) params = TIRTLtools::params
  df = .add_alleles(df) # add "*01" as allele for va and vb if necessary
  df = .filter_alleles(df, params = params) # remove alleles not found in the parameter data frame
  df = df %>% filter(nchar(cdr3a) > 5, nchar(cdr3b) > 5)
  #if(!"is_functional" %in% colnames(df))
  df = identify_non_functional_seqs(df)
  df = df %>% filter(is_functional) # remove seqs w/ stop codons (*) or frameshifts (_)
  df = as.data.frame(df) ## convert to standard data frame
  return(df)
}
