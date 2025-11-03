#' Write a tab-separated "thimble" text file for use with Stitchr
#'
#' `r lifecycle::badge('experimental')`
#'
#' @family stitchr
write_stitchr_tsv = function(df, output_name, output_folder = "") {
  cols = c("TCR_name", "TRAV", "TRAJ", "TRA_CDR3", "TRBV", "TRBJ", "TRB_CDR3",
           "TRAC", "TRBC", "TRA_leader", "TRB_leader", "Linker", "Link_order",
           "TRA_5_prime_seq", "TRA_3_prime_seq", "TRB_5_prime_seq", "TRB_3_prime_seq")
  if( sum(!cols %in% colnames(df)) > 0) {
    message("Thimble output requires a data frame with the following columns: TCR_name, TRAV, TRAJ, TRA_CDR3, TRBV, TRBJ, TRB_CDR3, TRAC, TRBC, TRA_leader, TRB_leader, Linker, Link_order, TRA_5_prime_seq, TRA_3_prime_seq, TRB_5_prime_seq, TRB_3_prime_seq.")
    message("Please use the paired_to_thimble_df() function to create a properly formatted data frame.")
    stop("Incorrectly formatted data frame.")
  }
  chk1 = sum( grepl("^m", df$TRAC) | grepl("^m", df$TRBC) ) > 0 ## check for mouse genes in TRAC and TRBC
  if(chk1) message("Note: You may need to run stitchr with the '-xg' option for 'additional/custom genes' since your TRAC or TRBC uses a mouse gene.")
  filename = paste(output_name, ".tsv", sep = "")
  ff = file.path(output_folder, filename)
  write.table(df, ff, row.names = FALSE, sep = "\t")
  return(invisible(NULL))
}
