rename_columns <- function(df, format = c("10X", "ParseBio"), rename_df = get_names_df_single_chain(), verbose = TRUE) {
  format = format[1]
  col_name = paste("name_", format, sep = "")
  # Validate inputs
  if (!is.data.frame(df)) stop("'df' must be a data frame.")
  if (!is.data.frame(rename_df)) stop("'rename_df' must be a data frame.")

  # Find which old_names are missing from df
  missing_cols <- rename_df[[col_name]][!rename_df[[col_name]] %in% colnames(df)]

  if (length(missing_cols) > 0 && verbose) {
    message(
      "The following column(s) were not found in the data frame and will be skipped: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  valid_renames <- rename_df[rename_df[[col_name]] %in% colnames(df), ]

  # Perform renaming
  colnames(df)[match(valid_renames[[col_name]], colnames(df))] <- valid_renames$name

  return(df)
}

get_names_df_single_chain = function() {
  rename_df <- tribble(
    ~name,      ~name_10X,   ~name_ParseBio,
    "barcode",  "barcode",   "cell_barcode",
    "chain",    "chain",     "locus",
    "cdr3_aa",  "cdr3",      "cdr3_aa",
    "cdr3_nt",  "cdr3_nt",   "cdr3",
    "v_gene",   "v_gene",    "v_call",
    "j_gene",   "j_gene",    "j_call",
    "reads",    "reads",     "read_count",
    "umis",     "umis",      "transcript_count"
  )
  return(rename_df)
}
