
## input is a list with "row_data" (TCR metadata), "col_data" (well metadata), and "matrix" (sparse matrix of TCRs x wells)
## col_data needs a column "well" with e.g. "A13", etc.
plot_nreads_plate = function(list) {
  df = list$col_data
  df$n_reads = Matrix::colSums(list$matrix)
  df$row = factor(substr(df$well, 0, 1), levels = rev(LETTERS[1:16]))
  df$column = factor(substr(df$well, 2,3), levels = as.character(1:24))
  ggplot(df) + geom_tile(aes(x=column, y=row, fill = n_reads)) +
    scale_y_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE) +
    scale_fill_gradient(low = "blue", high = "orange") +
    coord_fixed() +
    theme_classic()
}


plot_n_tcrs_one_well = function(list) {
  df = list$col_data
  n_wells_obs = Matrix::rowSums(list$matrix > 0)
  n_wells_one_obs = which(n_wells_obs == 1)
  mat_sub = list$matrix[n_wells_one_obs,]
  df$n_tcrs_one_well = Matrix::colSums(mat_sub > 0)
  df$row = factor(substr(df$well, 0, 1), levels = rev(LETTERS[1:16]))
  df$column = factor(substr(df$well, 2,3), levels = as.character(1:24))
  ggplot(df) + geom_tile(aes(x=column, y=row, fill = n_tcrs_one_well)) +
    scale_y_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE) +
    scale_fill_gradient(low = "blue", high = "orange") +
    coord_fixed() +
    theme_classic()
}


plot_n_tcrs = function(list) {
  df = list$col_data
  df$n_tcrs = Matrix::colSums(list$matrix > 0)
  df$row = factor(substr(df$well, 0, 1), levels = rev(LETTERS[1:16]))
  df$column = factor(substr(df$well, 2,3), levels = as.character(1:24))
  ggplot(df) + geom_tile(aes(x=column, y=row, fill = n_tcrs)) +
    scale_y_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE) +
    scale_fill_gradient(low = "blue", high = "orange") +
    coord_fixed() +
    theme_classic()
}


identify_bad_wells = function(alpha_list, beta_list, well_filter_thres = 0.5) {
  n_clones_alpha = Matrix::colSums(alpha_list$matrix > 0)
  n_clones_beta = Matrix::colSums(beta_list$matrix > 0)
  min_clones = round(well_filter_thres*mean(n_clones_alpha)) ## use alpha to match github pseudobulk code
  bad_wells_alpha = n_clones_alpha <= min_clones
  bad_wells_beta = n_clones_beta <= min_clones
  bad_wells = bad_wells_alpha | bad_wells_beta
  return(bad_wells)
}

identify_bad_wells_single = function(list, well_filter_thres = 0.5) {
  n_clones = Matrix::colSums(list$matrix > 0)
  min_clones = round(well_filter_thres*mean(n_clones)) ## thres*(avg num of clones found per well)
  bad_wells = n_clones <= min_clones
  print(paste("Bad wells:", paste(list$col_data$well[bad_wells], collapse = ", ")))
  return(bad_wells)
}

sparse_to_pseudobulk_single = function(mat, row_data) {
  if(nrow(mat) != nrow(row_data)) stop("Row data doesn't match size of matrix")
  keep_rows = sparseMatrixStats::rowSums2(mat) != 0
  mat = mat[keep_rows,]
  mat_without_zeros = mat
  #mat_without_zeros[mat_without_zeros == 0] = NA
  row_data = row_data[keep_rows,]
  rc = sparseMatrixStats::rowSums2(mat)
  rc_max = sparseMatrixStats::rowMaxs(mat)

  ### efficient code to compute medians, dropping zeros
  #rc_med = sparseMatrixStats::rowMedians(mat)
  M <- drop0(mat)
  split_x <- split(M@x, M@i + 1L) ## list of data for each row
  med <- rep(NA_real_, nrow(M))
  med[lengths(split_x) > 0] <- vapply(split_x, median, numeric(1))


  rc_med_w_zeros = sparseMatrixStats::rowMedians(mat)
  ncols = ncol(mat)
  nwells = apply(mat, 1, function(x) sum(x!=0) )
  df_out = row_data %>%
    as_tibble() %>%
    mutate(
      readCount = rc,
      readCount_max = rc_max,
      readCount_median = med,
      readCount_median_with_zeros = rc_med_w_zeros,
      n_wells = nwells,
      max_wells = ncols,
      readFraction = readCount/sum(readCount)
    ) %>%
    mutate(sem = sd(readFraction)/sqrt(ncols)) %>%
    arrange(desc(readFraction))
  return(df_out)
}

sparse_to_pseudobulk = function(alpha_list, beta_list) {
  good_wells = !identify_bad_wells(alpha_list, beta_list)
  mat_alpha = alpha_list$matrix[,good_wells]
  row_data_alpha = alpha_list$row_data
  mat_beta = beta_list$matrix[,good_wells]
  row_data_beta = beta_list$row_data

  pb_alpha = sparse_to_pseudobulk_single(mat_alpha, row_data_alpha)
  pb_beta = sparse_to_pseudobulk_single(mat_beta, row_data_beta)
  return(list(alpha = pb_alpha, beta = pb_beta))
}

