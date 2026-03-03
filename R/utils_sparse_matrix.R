# sparse matrix helpers


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

extract_col_dense <- function(mat, j) {
  col_vec <- numeric(nrow(mat))
  p_start <- mat@p[j] + 1
  p_end   <- mat@p[j + 1]
  if (p_end >= p_start) {
    col_vec[mat@i[p_start:p_end] + 1] <- mat@x[p_start:p_end]
  }
  col_vec
}

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
