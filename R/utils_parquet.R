# parquet helpers

lazy_load_parquet = function(file) {
  con <- duckdb::dbConnect(duckdb::duckdb())
  call = paste("read_parquet(", "'", file, "'", ")", sep = "")
  #tbl <- tbl(con, "read_parquet('alpha_test.parquet')")
  df <- dplyr::tbl(con, call)
  return(df)
}

get_nrows_parquet = function(file) {
  con <- duckdb::dbConnect(duckdb::duckdb())
  nrows <- DBI::dbGetQuery(
    con,
    sprintf("SELECT COUNT(*) FROM read_parquet('%s')", file)
  )[[1]]
  return(nrows)
}

extract_rows_by_index_parquet <- function(parquet_path, row_indices) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  # Register index vector as a temp table — much faster than huge IN() clause
  idx_df <- data.frame(rn = as.integer(row_indices) - 1L)
  DBI::dbWriteTable(con, "row_idx", idx_df, temporary = TRUE)

  query <- sprintf("
    SELECT * EXCLUDE (rn)
    FROM (
      SELECT *, (row_number() OVER () - 1)::INTEGER AS rn
      FROM read_parquet('%s')
    ) t
    INNER JOIN row_idx USING (rn)
  ", parquet_path)

  DBI::dbGetQuery(con, query)
}

get_row_parquet = function(file, col, value, max_rows = 1) {
  con <- duckdb::dbConnect(duckdb::duckdb())
  #specific_rows <- dbGetQuery(con, "SELECT * FROM read_parquet('file.parquet') WHERE column_name = 'value'")
  #call = paste("read_parquet(", "'", file, "'", ")", sep = "")
  call = paste("SELECT * FROM read_parquet('", file, "') WHERE ", col ," = '", value, "'", sep = "")
  df = DBI::dbGetQuery(con,call)
  nrows = min(max_rows, nrow(df))
  #df2 <- dplyr::tbl(con, call)
  return(df[1:nrows,])
}

get_row_data = function(se, backend = c("duckdb", "arrow", "duckplyr")) {
  backend = backend[1]
  file = metadata(se)$row_data_file
  if(backend == "arrow") df = arrow::open_dataset(file)
  if(backend == "duckdb") {
    con <- duckdb::dbConnect(duckdb::duckdb())
    call = paste("read_parquet(", "'", file, "'", ")", sep = "")
    #tbl <- tbl(con, "read_parquet('alpha_test.parquet')")
    df <- dplyr::tbl(con, call)
  }
  if(backend == "duckplyr") df = duckplyr::read_parquet_duckdb(file)
  return(df)
}

find_row_index <- function(parquet_path, col, value) {
  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))

  query <- sprintf("
    SELECT rn + 1 AS row_index  -- convert back to 1-based
    FROM (
      SELECT row_number() OVER () - 1 AS rn, %s
      FROM read_parquet('%s')
    )
    WHERE %s = ?
  ", col, parquet_path, col)

  result <- DBI::dbGetQuery(con, query, params = list(value))
  result$row_index
}
