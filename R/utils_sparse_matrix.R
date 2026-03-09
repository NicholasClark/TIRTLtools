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
#' @noRd
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
#' @noRd
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

## for correlations of all vs. all
sparse_col_cor_all <- function(X, Y) {
  # X and Y must be dgCMatrix with the same number of rows
  stopifnot(nrow(X) == nrow(Y))

  n <- nrow(X)

  # --- Efficient sparse column stats ---
  col_means <- function(M) Matrix::colMeans(M)

  col_sd <- function(M) {
    # Var = E[x^2] - E[x]^2, computed sparsely
    means  <- col_means(M)
    means2 <- Matrix::colMeans(M^2)
    sqrt(pmax(means2 - means^2, 0) * (n / (n - 1)))  # unbiased (Bessel's correction)
  }

  # --- Center columns (returns a regular dense matrix only for cross-product) ---
  # We compute cov(X, Y) = (1/(n-1)) * t(X_centered) %*% Y_centered sparsely

  mx <- col_means(X)
  my <- col_means(Y)
  sx <- col_sd(X)
  sy <- col_sd(Y)

  # Efficient cross-product using sparse arithmetic:
  # cov(xi, yj) = [sum(xi * yj) - n * mean(xi) * mean(yj)] / (n - 1)
  # crossprod gives t(X) %*% Y  (p x q), sum of elementwise products per column pair
  XtY <- Matrix::crossprod(X, Y)  # p x q sparse/dense matrix

  # Outer product of means scaled by n
  mean_correction <- outer(mx, my) * n  # p x q dense

  cov_mat <- (as.matrix(XtY) - mean_correction) / (n - 1)  # p x q

  # Correlation = cov / (sd_x * sd_y)
  cor_mat <- cov_mat / outer(sx, sy)  # p x q

  # Clip to [-1, 1] to avoid numerical issues
  cor_mat <- pmin(pmax(cor_mat, -1), 1)

  # --- t-statistic and p-value ---
  # t = r * sqrt(n - 2) / sqrt(1 - r^2)
  t_mat <- cor_mat * sqrt(n - 2) / sqrt(1 - cor_mat^2)

  # Two-sided p-value from t-distribution with df = n - 2
  p_mat <- 2 * pt(-abs(t_mat), df = n - 2)

  list(
    correlation = cor_mat,
    t_statistic = t_mat,
    p_value     = p_mat
  )
}
