#' Convert paired TCRs to a "thimble" data frame for stitching together TCR nucleotide coding sequences with the Stitchr Python package
#'
#' @family stitchr
#'
create_thimble_df = function(df, preset = c("default", "none"), TCR_names = "TCR",
                             exclude_non_functional = TRUE, remove_duplicates = TRUE, verbose = TRUE,
                             Linker = NULL, Link_order = NULL, TRAC = NULL, TRBC = NULL,
                             TRA_leader = NULL, TRB_leader = NULL, TRA_5_prime_seq = NULL,
                             TRA_3_prime_seq = NULL, TRB_5_prime_seq = NULL, TRB_3_prime_seq = NULL) {
  preset = preset[1]
  if(remove_duplicates) {
    dupes = duplicated(paste(df$va, df$vb, df$ja, df$jb, df$cdr3a, df$cdr3b))
    df = df[!dupes,]
    msg0 = paste("Note: Removed",sum(dupes), "duplicate TCRs.")
    if(verbose) message(msg0)
  }
  if(exclude_non_functional) {
    df = identify_non_functional_seqs(df)
    sum_nonfunc = sum(!df$is_functional)
    msg1 = paste("Note: Removed", sum_nonfunc, "TCRs with non-functional CDR3 sequences.")
    if(verbose) message(msg1)
    df = df %>% filter(is_functional)
  }
  cols = c("TCR_name", "TRAV", "TRAJ", "TRA_CDR3", "TRBV", "TRBJ", "TRB_CDR3",
           "TRAC", "TRBC", "TRA_leader", "TRB_leader", "Linker", "Link_order",
           "TRA_5_prime_seq", "TRA_3_prime_seq", "TRB_5_prime_seq", "TRB_3_prime_seq")
  df_tmp = data.frame(matrix("", ncol = length(cols), nrow = nrow(df)))
  colnames(df_tmp) = cols
  df_tmp$TRAV = df$va
  df_tmp$TRAJ = df$ja
  df_tmp$TRA_CDR3 = df$cdr3a
  df_tmp$TRBV = df$vb
  df_tmp$TRBJ = df$jb
  df_tmp$TRB_CDR3 = df$cdr3b
  df_tmp$TCR_name = paste("TCR", 1:nrow(df_tmp), sep = "_")
  if(preset == "default") {
    df_tmp$Linker = "T2A"
    df_tmp$Link_order = "BA"
    df_tmp$TRAC = "mTRAC*01"
    df_tmp$TRBC = "mTRBC2*01"
  }
  ### add in code to set other parameters
  input_cols = c("TRAC", "TRBC", "TRA_leader", "TRB_leader", "Linker", "Link_order",
                 "TRA_5_prime_seq", "TRA_3_prime_seq", "TRB_5_prime_seq", "TRB_3_prime_seq")
  for(col in input_cols) {
    if(!is.null(get(col))) df_tmp[[col]] = get(col)
  }
  if(length(TCR_names) == 1) TCR_names = paste(TCR_names, 1:(dim(df_tmp)[1]), sep = "_")
  df_tmp$TCR_name = TCR_names
  df_tmp$TCR_name = make.unique(df_tmp$TCR_name, sep = "_")
  chk1 = sum( grepl("^m", df_tmp$TRAC) | grepl("^m", df_tmp$TRBC) ) > 0 ## check for mouse genes in TRAC and TRBC
  if(chk1) message("Note: You may need to run stitchr with the '-xg' option for 'additional/custom genes' since your TRAC or TRBC uses a mouse gene.")
  return(df_tmp)
}
