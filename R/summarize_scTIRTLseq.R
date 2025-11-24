# Summarize single-cell TIRTLseq data
#
# @description
# `r lifecycle::badge('experimental')`
# Returns a summary of single-cell TIRTLseq data, printing the total number of
# wells with data, number of wells missing an alpha/beta chain, number of wells
# with a second alpha chain, etc.
#
# @param df a data frame with single-cell TIRTLseq data
#
#
# @family single-cell
summarize_scTIRTLseq = function(df) {
  n_wells_total = dim(df)[1]
  n_wells_both = sum( (!is.na(df$aaSeqCDR3u_alpha)) & (!is.na(df$aaSeqCDR3u_beta)) )
  percent_both = 100*(n_wells_both/n_wells_total) %>% round(3)
  n_wells_no_alpha = sum(is.na(df$aaSeqCDR3u_alpha))
  n_wells_no_beta = sum(is.na(df$aaSeqCDR3u_beta))
  percent_no_alpha = 100*(n_wells_no_alpha/n_wells_total) %>% round(3)
  percent_no_beta = 100*(n_wells_no_beta/n_wells_total) %>% round(3)

  n_wells_second_alpha = sum(!is.na(df$aaSeqCDR3u_alpha_second))
  percent_second_alpha = 100*(n_wells_second_alpha/n_wells_total) %>% round(3)

  msg1 = paste(n_wells_total, " wells of data", sep = "")
  msg2 = paste(n_wells_both, " wells (", percent_both, "%) with both alpha and beta chains", sep = "")
  msg3 = paste(n_wells_no_alpha, " wells (", percent_no_alpha, "%) missing an alpha chain", sep = "")
  msg4 = paste(n_wells_no_beta, " wells (", percent_no_beta, "%) missing an alpha chain", sep = "")
  msg5 = paste(n_wells_second_alpha, " wells (", percent_second_alpha, "%) with second alpha chain")

  message(msg1)
  message(msg2)
  message(msg3)
  message(msg4)
  message(msg5)
  return(invisible(NULL))
}
