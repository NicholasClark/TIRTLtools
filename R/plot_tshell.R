plot_tshell_test = function(
    well_data,
    nuc_seq,
    chain,
    n_plot = 9,
    #plot_overlap_only = TRUE,
    interactive = FALSE
) {
  if(!chain %in% c("alpha", "beta")) stop("'chain' needs to be 'alpha' or 'beta'")


  wa_all = Matrix::colSums(well_data$alpha != 0)
  wb_all = Matrix::colSums(well_data$beta != 0)

  alpha_mat = well_data$alpha
  beta_mat = well_data$beta

  alpha_mat_rf = normalize_rows(alpha_mat)
  beta_mat_rf = normalize_rows(beta_mat)

  is_lazy = !"data.frame" %in% class(well_data$alpha_meta)
  if(!is_lazy) {
    if(chain == "alpha") idx = which(ll_alpha$row_data$targetSequences == nuc_seq)[1]
    if(chain == "beta") idx = which(ll_beta$row_data$targetSequences == nuc_seq)[1]
  } else {
    if(chain == "alpha") idx = find_row_index(well_data$alpha_meta_file, "targetSequences", nuc_seq)[1]
    if(chain == "beta") idx = find_row_index(well_data$beta_meta_file, "targetSequences", nuc_seq)[1]
  }

  #cor_df = ll_beta$row_data

  wa1 = wa_all[idx]
  alpha_meta = extract_rows_by_index_parquet(well_data$alpha_meta_file, idx)
  alpha_vec = extract_col_dense(alpha_mat_rf, idx)

  # cor_df = well_data[[paste("rows", other, sep = "_")]]
  # mat_dense_t = t(as.matrix(mat_other))
  # well_sums = Matrix::rowSums(mat_dense_t)
  # mat_dense_t = mat_dense_t/well_sums
  #
  # cor_df$w_other = Matrix::rowSums(mat_other != 0)
  # cor_df$w_chain = sum(!is.na(df_single$readFraction))

  #mat_dense_t[mat_dense_t <= 0] = NA
  #mat_log10 = log10(dense_mat)

  cor_df = sparse_col_cor(beta_mat_rf, alpha_vec)
  #cor_df = bind_cols(cor_df, res)

  #cor_df$r_log_log = cor(mat_log10, log10(df_single$readFraction), method = "pearson", use = "pairwise.complete.obs")
  cor_df$wij = sparse_overlap(beta_mat_rf, alpha_vec)
  cor_df$wa = wa1
  cor_df$wb = wb_all
  cor_df$p_adj = p.adjust(cor_df$p, method = "BH")

  keep_indices = which(cor_df$wij > 0.5*max(cor_df$wij))
  cor_df2 = cor_df[keep_indices,]

  beta_meta2 = extract_rows_by_index_parquet(well_data$beta_meta_file, row_indices = keep_indices)
  cor_df2 = bind_cols(beta_meta2, cor_df2) %>% as_tibble() %>%
    select(chain, row_id, aaSeqCDR3, r, t, p, p_adj, wa, wb, wij, n, targetSequences, v,j, everything()) %>%
    arrange(desc(r))


  cor_df3 = cor_df2[1:n_plot,]

  meta_gg = extract_rows_by_index_parquet(well_data$beta_meta_file, cor_df3$row_id) %>%
    mutate(chain_name = paste("beta", row_id))
  mat_gg = as.matrix(beta_mat_rf[,cor_df3$row_id])

  colnames(mat_gg) = meta_gg$targetSequences
  rownames(mat_gg) = well_data$col_meta$well
  gg_df = lapply(1:n_plot, function(i) {
    df_tmp = tibble(readFraction1 = mat_gg[,i],
                    readFraction2 = alpha_vec,
                    well = well_data$col_meta$well,
                    chain_name = meta_gg$chain_name[i]
    )
  }) %>% bind_rows()
  gg_df$chain_name = factor(gg_df$chain_name, levels = meta_gg$chain_name)
  gg_df$readFraction1[gg_df$readFraction1 == 0] = NA
  #print(cor_df3)
  gg = ggplot(gg_df, aes(x=readFraction1, y=readFraction2, label = well)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~chain_name, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  #print(gg)
  gg2 = ggplot(cor_df2, aes(x=r, y=-log10(p), aa = aaSeqCDR3, v = v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  gg3 = ggplot(cor_df2, aes(x=r, y=-log10(p_adj), aa = aaSeqCDR3, v= v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  if(interactive) {
    gg = plotly::ggplotly(gg)
    gg2 = plotly::ggplotly(gg2)
    gg3 = plotly::ggplotly(gg3)
  }

  return(list(cor_df = cor_df, cor_df2 = cor_df2, cor_df3 = cor_df3, gg = gg, gg2 = gg2, gg3 = gg3))
}

plot_tshell = function(
    ll_alpha,
    ll_beta,
    nuc_seq,
    chain,
    n_plot = 10,
    #plot_overlap_only = TRUE,
    interactive = FALSE
) {
  if(!chain %in% c("alpha", "beta")) stop("'chain' needs to be 'alpha' or 'beta'")

  ll_alpha$row_data$wa = Matrix::colSums(ll_alpha$matrix!= 0)
  ll_beta$row_data$wb = Matrix::colSums(ll_beta$matrix!= 0)

  alpha_mat = ll_alpha$matrix
  beta_mat = ll_beta$matrix

  alpha_mat_rf = normalize_rows(alpha_mat)
  beta_mat_rf = normalize_rows(beta_mat)

  #if("data.frame" %in% class(test$rows_alpha))
  if(chain == "alpha") idx = which(ll_alpha$row_data$targetSequences == nuc_seq)[1]
  if(chain == "beta") idx = which(ll_beta$row_data$targetSequences == nuc_seq)[1]

  cor_df = ll_beta$row_data

  idx = which(ll_alpha$row_data$targetSequences == nuc_seq)[1]
  wa1 = ll_alpha$row_data$wa[idx]
  alpha_vec = extract_col_dense(alpha_mat_rf, idx)

  # cor_df = well_data[[paste("rows", other, sep = "_")]]
  # mat_dense_t = t(as.matrix(mat_other))
  # well_sums = Matrix::rowSums(mat_dense_t)
  # mat_dense_t = mat_dense_t/well_sums
  #
  # cor_df$w_other = Matrix::rowSums(mat_other != 0)
  # cor_df$w_chain = sum(!is.na(df_single$readFraction))

  #mat_dense_t[mat_dense_t <= 0] = NA
  #mat_log10 = log10(dense_mat)

  res = sparse_col_cor(beta_mat_rf, alpha_vec)
  cor_df = bind_cols(cor_df, res)

  #cor_df$r_log_log = cor(mat_log10, log10(df_single$readFraction), method = "pearson", use = "pairwise.complete.obs")
  cor_df$wij = sparse_overlap(beta_mat_rf, alpha_vec)
  cor_df$wa = wa1
  cor_df$p_adj = p.adjust(cor_df$p, method = "BH")

  cor_df = cor_df %>% select(chain, row_id, aaSeqCDR3, r, t, p, p_adj, wa, wb, wij, n, targetSequences, v,j, everything())

  cor_df2 = cor_df %>% filter(wij > 0.5*max(cor_df$wij)) %>% arrange(desc(r))

  cor_df3 = cor_df2[1:n_plot,]

  meta_gg = ll_beta$row_data[match(cor_df3$row_id, ll_beta$row_data$row_id),] %>%
    mutate(chain_name = paste("beta", row_id))
  mat_gg = as.matrix(beta_mat_rf[,cor_df3$row_id])

  colnames(mat_gg) = meta_gg$targetSequences
  rownames(mat_gg) = ll_beta$col_data$well
  gg_df = lapply(1:n_plot, function(i) {
    df_tmp = tibble(readFraction1 = mat_gg[,i],
                    readFraction2 = alpha_vec,
                    well = ll_beta$col_data$well,
                    chain_name = meta_gg$chain_name[i]
    )
  }) %>% bind_rows()
  gg_df$chain_name = factor(gg_df$chain_name, levels = meta_gg$chain_name)
  gg_df$readFraction1[gg_df$readFraction1 == 0] = NA
  #print(cor_df3)
  gg = ggplot(gg_df, aes(x=readFraction1, y=readFraction2, label = well)) +
    geom_point(alpha = 0.3) +
    facet_wrap(~chain_name, scales = "free") +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw()
  #print(gg)
  gg2 = ggplot(cor_df2, aes(x=r, y=-log10(p), aa = aaSeqCDR3, v = v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  gg3 = ggplot(cor_df2, aes(x=r, y=-log10(p_adj), aa = aaSeqCDR3, v= v, j = j, wa = wa, wb = wb, wij = wij)) + geom_point(alpha = 0.5) + theme_bw()
  if(interactive) {
    gg = plotly::ggplotly(gg)
    gg2 = plotly::ggplotly(gg2)
    gg3 = plotly::ggplotly(gg3)
  }

  return(list(cor_df = cor_df, cor_df2 = cor_df2, cor_df3 = cor_df3, gg = gg, gg2 = gg2, gg3 = gg3))
}

extract_col_dense <- function(mat, j) {
  col_vec <- numeric(nrow(mat))
  p_start <- mat@p[j] + 1
  p_end   <- mat@p[j + 1]
  if (p_end >= p_start) {
    col_vec[mat@i[p_start:p_end] + 1] <- mat@x[p_start:p_end]
  }
  col_vec
}

get_nrows_parquet = function(file) {
  con <- duckdb::dbConnect(duckdb::duckdb())
  nrows <- DBI::dbGetQuery(
    con,
    sprintf("SELECT COUNT(*) FROM read_parquet('%s')", file)
  )[[1]]
  return(nrows)
}

# nuc_seq = "TGTGTCGGGATGGATAGCAGCTATAAATTGATCTTC" ## 8 partners
# nuc_seq = "TGTGCATCTCCGCCTCTTAACAAATTTTACTTT" # only 3 wells... R near 1, problem?
# nuc_seq = "TGTGCTGTGAGTGACAAAAGCGACTACAAGCTCAGCTTT" # very clear good one
# plot_tshell(ll_alpha, ll_beta, nuc_seq="TGTGCTGTGAGTGACAAAAGCGACTACAAGCTCAGCTTT", chain = "alpha")
# plot_tshell(ll_alpha, ll_beta, nuc_seq="TGTGTCGGGATGGATAGCAGCTATAAATTGATCTTC", chain = "alpha")
# plot_tshell(ll_alpha, ll_beta, nuc_seq="TGTGCATCTCCGCCTCTTAACAAATTTTACTTT", chain = "alpha")
#well_data = load_wells_counts_rds("/Volumes/thomagrp/home/common/nclark52/temp/MVP175_Bleakley_miH/TCR_clones_UD18_rds", prefix = "MVP175", lazy = TRUE)




# sparse_col_cor_old <- function(X, y) {
#   stopifnot(is(X, "sparseMatrix"), length(y) == nrow(X))
#
#   n <- nrow(X)
#
#   # --- y stats ---
#   y_c  <- y - mean(y)          # centered y; sum(y_c) = 0
#   sd_y <- sqrt(sum(y_c^2) / (n - 1))
#
#   # --- Numerator: t(X_centered) %*% y_c ---
#   # Since sum(y_c) = 0:
#   # sum_i (X_ij - mean_j)(y_i - mean_y) = sum_i X_ij * y_c_i  (cross term vanishes)
#   num <- as.numeric(Matrix::crossprod(X, y_c))   # avoids forming t(X) densely
#
#   # --- Column std of X (sparse-aware, O(nnz)) ---
#   col_means  <- Matrix::colMeans(X)
#   col_sq_means <- Matrix::colMeans(X ^ 2)        # X^2 stays sparse for sparse X
#   sd_x <- sqrt(pmax(col_sq_means - col_means^2, 0) * n / (n - 1))
#
#   # --- Pearson r ---
#   r <- num / ((n - 1) * sd_x * sd_y)
#   r[sd_x == 0] <- NA            # constant columns → NA
#
#   setNames(r, colnames(X))
#   return(r)
# }

#' Correlate a dense vector y with every column of a sparse matrix X,
#' returning r, t-statistic, and p-value
#' @param X  dgCMatrix of shape (n x p)
#' @param y  numeric vector of length n
#' @return   data.frame with columns: r, t, p, n_overlap
sparse_col_cor <- function(X, y, use_overlap = FALSE) {
  stopifnot(is(X, "sparseMatrix"), length(y) == nrow(X))

  n <- nrow(X)

  # --- y stats ---
  y_c  <- y - mean(y)
  sd_y <- sqrt(sum(y_c^2) / (n - 1))

  # --- Numerator: avoids centering X (cross-term vanishes since sum(y_c)=0) ---
  num <- as.numeric(Matrix::crossprod(X, y_c))

  # --- Column std of X ---
  col_means    <- Matrix::colMeans(X)
  col_sq_means <- Matrix::colMeans(X^2)
  sd_x <- sqrt(pmax(col_sq_means - col_means^2, 0) * n / (n - 1))

  # --- Pearson r ---
  r <- num / ((n - 1) * sd_x * sd_y)
  r[sd_x == 0] <- NA

  # --- Effective n (optionally use overlap of nonzeros) ---
  if (use_overlap) {
    # Reuse sparse overlap logic: how many rows are nonzero in both X[:,j] and y
    y_sp    <- as(as(y, "sparseVector"), "CsparseMatrix")
    y_sp@x  <- rep(1, length(y_sp@x))
    X_bin   <- X; X_bin@x[] <- 1
    n_eff   <- as.integer(Matrix::crossprod(y_sp, X_bin))
  } else {
    n_eff <- rep(n, ncol(X))
  }

  # --- t-statistic: r * sqrt((n-2) / (1-r^2)) ---
  t_stat <- r * sqrt((n_eff - 2) / pmax(1 - r^2, .Machine$double.eps))

  # --- Two-tailed p-value from t distribution ---
  p_val <- 2 * pt(abs(t_stat), df = n_eff - 2, lower.tail = FALSE)

  data.frame(
    r       = r,
    t       = t_stat,
    p       = p_val,
    n       = n_eff,
    row.names = colnames(X)
  )
}


#' Normalize rows of a dgCMatrix by their row sums
#' @param X  dgCMatrix
#' @return   dgCMatrix with rows summing to 1 (zero-sum rows → NA or 0)
normalize_rows <- function(X) {
  rs <- Matrix::rowSums(X)
  if (any(rs == 0)) stop("Zero-sum rows found: ", paste(which(rs == 0), collapse = ", "))
  # Efficient: multiply by diagonal matrix from the left
  # Equivalent to sweep(X, 1, rs, "/") but stays sparse
  X_norm <- Diagonal(x = 1 / rs) %*% X
  return(X_norm)
}

sparse_overlap <- function(X, y) {
  stopifnot(is(X, "sparseMatrix"), length(y) == nrow(X))

  # Convert y to a sparse column vector — only nonzero entries matter
  y_sparse <- as(as(y, "sparseVector"), "CsparseMatrix")  # n x 1 dgCMatrix

  # Binarize both: nonzero -> 1
  # This replaces stored values without changing sparsity structure
  X@x[]       <- rep(1, length(X@x))
  y_sparse@x  <- rep(1, length(y_sparse@x))

  # crossprod: t(y_sparse) %*% X  gives a 1xp matrix of overlap counts
  # Only nnz(y) rows of X are ever touched
  as.integer(Matrix::crossprod(y_sparse, X))
}

# extract_rows_by_index_parquet <- function(parquet_path, row_indices) {
#   # DuckDB uses 1-based row numbers via rowid (actually 0-based internally)
#   # We'll use a generated row number approach
#
#   con <- DBI::dbConnect(duckdb::duckdb())
#   on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
#
#   # Convert indices to a comma-separated string for SQL IN clause
#   indices_str <- paste(row_indices - 1L, collapse = ", ")  # convert to 0-based
#
#   query <- sprintf("
#     SELECT * EXCLUDE (rn)
#     FROM (
#       SELECT *, row_number() OVER () - 1 AS rn
#       FROM read_parquet('%s')
#     )
#     WHERE rn IN (%s)
#   ", parquet_path, indices_str)
#
#   DBI::dbGetQuery(con, query)
# }

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

