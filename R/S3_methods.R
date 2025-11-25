# Create a TIRTLseqDataSet object (S3)
# @export
# TIRTLseqDataSet <- function(meta, data, flags = list()) {
#   stopifnot(
#     is.data.frame(meta),
#     is.list(data),
#     !is.null(names(data))
#   )
#
#   # Optional: sanity checks on per-sample structure
#   ok <- vapply(data, function(x) {
#     is.list(x) &&
#       all(c("alpha", "beta", "paired") %in% names(x)) &&
#       all(vapply(x[c("alpha", "beta", "paired")], is.data.frame, logical(1)))
#   }, logical(1))
#
#   if (!all(ok)) {
#     stop("Each element of 'data' must be a list with data.frames 'alpha', 'beta', and 'paired'.")
#   }
#
#   structure(
#     list(
#       meta  = meta,
#       data  = data,
#       flags = flags
#     ),
#     class = "TIRTLobj"
#   )
# }

#' @export
print.TIRTLseqDataSet <- function(x, ...) {
  n_samples <- length(x$data)
  cat("<TIRTLseqDataSet>\n")
  cat("  Number of samples:", n_samples, "\n")
  if (n_samples == 0L) {
    cat("  (no samples)\n")
    return(invisible(x))
  } else {
    sample_names <- names(x$data)
    cat("  Samples:", sample_names)
  }
  invisible(x)
}


summary.TIRTLseqDataSet = function(x, ...) {
  n_samples <- length(x$data)
  sample_names <- names(x$data)
  cat("<TIRTLseqDataSet>\n")
  cat("  Number of samples:", n_samples, "\n")

  # total TCR counts
  alpha_rows  <- integer(n_samples)
  beta_rows   <- integer(n_samples)
  paired_rows <- integer(n_samples)

  # unique TCR counts
  paired_rows2 <- integer(n_samples)

  for (i in seq_len(n_samples)) {
    s <- x$data[[i]]
    alpha_rows[i]  <- nrow(s$alpha)
    beta_rows[i]   <- nrow(s$beta)
    paired_rows[i] <- nrow(s$paired)

    s_unique = remove_duplicates(s)
    paired_rows2[i] <- nrow(s_unique$paired)
  }
  cat("  Per-sample TCR counts (alpha / beta / paired):\n")
  summary_df <- data.frame(
    sample_id      = sample_names,
    `Paired TCRs` = paired_rows,
    `Unique Paired TCRs` = paired_rows2,
    `Alpha Chains`  = alpha_rows,
    `Beta Chains`   = beta_rows,
    row.names   = NULL,
    check.names = FALSE
  )

  print(tibble(summary_df))

  # max_show <- 10L
  # cat("  Per-sample row counts (alpha / beta / paired):\n")
  # if (n_samples > max_show) {
  #   print(utils::head(summary_df, max_show), row.names = FALSE)
  #   cat("  ...", n_samples - max_show, "more samples\n")
  # } else {
  #   print(summary_df, row.names = FALSE)
  # }
}
