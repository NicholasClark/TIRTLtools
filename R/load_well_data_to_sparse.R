
# load_well_data_to_sparse = function(folder_path, wells = get_well_subset(1:16,1:24),
#                           well_pos=3, chain = c("both", "alpha", "beta"), nproc = 1L,
#                           columns = NULL,
#                           max_files = Inf) {
#   chain = chain[1]
#
#   files = list.files(path = folder_path, full.names = F)
#   files_full = file.path(folder_path, files)
#   files_no_ending = gsub("\\.tsv", "", files)
#   #files_no_ending = gsub("\\.", "_", files_no_ending) ## replace periods with underscores if necessary
#
#   wells_all = sapply(strsplit(files_no_ending, "_"), function(x) x[[well_pos]])
#   ## make sure wells in the string are 2-3 characters long and only alphanumeric -- necessary for some data
#   wells_all = sapply(wells_all, function(x) substr(x, 1, 3))
#   wells_all = gsub("[^[:alnum:] ]", "", wells_all)
#
#   chains_all = sapply(strsplit(files_no_ending, "_"), function(x) rev(x)[[1]])
#   chains_all[chains_all=="TRA"] = "alpha"
#   chains_all[chains_all=="TRB"] = "beta"
#
#   wells_keep = wells_all %in% wells
#   if(chain %in% c("alpha","beta")) {
#     chains_keep = chains_all == chain
#   } else {
#     chains_keep = TRUE
#   }
#   keep = which(wells_keep & chains_keep)
#
#   meta_tmp = tibble(well = wells_all[keep], chain = chains_all[keep], file_short = files[keep], file_full = files_full[keep])
#   print(paste(dim(meta_tmp)[1], "files total"))
#   n_files = min(max_files, dim(meta_tmp)[1])
#   print(paste("Loading", n_files, "files"))
#   ll = parallel::mclapply(1:n_files, function(i) {
#     print(i)
#     file_full_tmp = meta_tmp$file_full[i]
#     file_short_tmp = meta_tmp$file_short[i]
#     well_tmp = meta_tmp$well[i]
#     chain_tmp = meta_tmp$chain[i]
#     #print(file_full_tmp)
#     df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
#       mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
#       mutate(well_chain = paste(well, chain))
#     return(df_tmp)
#   }, mc.cores = as.integer(nproc))
#   #names(ll) = files_no_ending[keep]
#   return(ll)
# }


load_well_data_to_sparse_multi = function(
    folder_path_list,
    sample_names_list,
    wells_list,
    well_pos_list,
    chain = c("alpha", "beta"),
    nproc = data.table::getDTthreads(),
    columns = NULL,
    max_files = Inf) {
  chain = chain[1]
  ll = lapply(1:length(folder_path_list), function(i) {
    print(i)
    list_tmp = load_well_data(folder_path = folder_path_list[[i]], wells = wells_list[[i]], well_pos = well_pos_list[[i]],
                   chain = chain, nproc = nproc, columns = columns, max_files = max_files) %>%
      bind_rows() %>% mutate(sample_id = sample_names_list[[i]], sample_num = i)
  })
  out = well_data_df_to_sparse_multi(ll)
  return(out)
}

load_well_data_to_sparse_single = function(
    folder_path,
    wells = get_well_subset(1:16,1:24),
    well_pos=3,
    chain = c("alpha", "beta"),
    file_type = c(".tsv", ".parquet"),
    parallel = FALSE,
    nproc = data.table::getDTthreads(),
    columns = NULL,
    max_files = Inf,
    periods_to_underscores = TRUE
    ) {
  tictoc::tic()
  args <- as.list(environment())
  df = do.call("load_well_data_to_df_single", args)
  ll = well_data_df_to_sparse(df, cols = wells)
  tictoc::toc()
  return(ll)

}

load_well_data_to_df_single = function(
    folder_path,
    wells = get_well_subset(1:16,1:24),
    well_pos=3,
    chain = c("alpha", "beta"),
    file_type = c(".tsv", ".parquet"),
    parallel = FALSE,
    nproc = data.table::getDTthreads(),
    columns = NULL,
    max_files = Inf,
    periods_to_underscores = TRUE
    ) {

  file_type = file_type[1]
  chain = chain[1]
  if(!chain %in% c("alpha", "beta")) stop("'chain' must be 'alpha' or 'beta'.")
  if(file_type == ".tsv") {
    chk1 = list.files(path = folder_path, pattern = "\\.gz$")
    if(length(chk1) > 0) warning("Folder contains .gz files -- check to make sure they have all been extracted to .tsv")
    if(chain == "alpha") files = list.files(path = folder_path, full.names = F, pattern = "_TRA\\.tsv$")
    if(chain == "beta") files = list.files(path = folder_path, full.names = F, pattern = "_TRB\\.tsv$")
    files_full = file.path(folder_path, files)
    files_no_ending = gsub("\\.tsv", "", files)
  } else if(file_type == ".parquet") {
    if(chain == "alpha") files = list.files(path = folder_path, full.names = F, pattern = "_TRA\\.parquet$")
    if(chain == "beta") files = list.files(path = folder_path, full.names = F, pattern = "_TRB\\.parquet$")
    files_full = file.path(folder_path, files)
    files_no_ending = gsub("\\.parquet", "", files)
  }

  if(periods_to_underscores) files_no_ending = gsub("\\.", "_", files_no_ending) ## replace periods with underscores if necessary

  wells_all = sapply(strsplit(files_no_ending, "_"), function(x) x[[well_pos]])
  ## make sure wells in the string are 2-3 characters long and only alphanumeric -- necessary for some data
  wells_all = sapply(wells_all, function(x) substr(x, 1, 3))
  wells_all = gsub("[^[:alnum:] ]", "", wells_all)

  wells_keep = wells_all %in% wells
  keep = which(wells_keep)

  meta_tmp = tibble(well = wells_all[keep], chain = chain, file_short = files[keep], file_full = files_full[keep])
  print(paste(dim(meta_tmp)[1], "files total"))
  n_files = min(max_files, dim(meta_tmp)[1])
  print(paste("Loading", n_files, "files"))
  cat("\n")
  df = my_apply(1:n_files, function(i) {
    cat("|")
    file_full_tmp = meta_tmp$file_full[i]
    file_short_tmp = meta_tmp$file_short[i]
    well_tmp = meta_tmp$well[i]
    chain_tmp = meta_tmp$chain[i]
    #print(file_full_tmp)
    if(file_type == ".tsv") {
      df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
        mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
        mutate(well_chain = paste(well, chain))
    } else if(file_type == ".parquet") {
      df_tmp = arrow::read_parquet(file_full_tmp) %>%
        mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
        mutate(well_chain = paste(well, chain))
      if(!is.null(columns)) df_tmp = df_tmp[,columns, drop = FALSE]
    }
    return(df_tmp)
  }, parallel = parallel, mc.cores = nproc) %>% bind_rows()  %>%
    dplyr::rename(j = allJHitsWithScore, v = allVHitsWithScore) %>%
    mutate(j = gsub("\\(.*", "", j), v = gsub("\\(.*", "", v))

  # df = parallel::mclapply(1:n_files, function(i) {
  #   cat("|")
  #   file_full_tmp = meta_tmp$file_full[i]
  #   file_short_tmp = meta_tmp$file_short[i]
  #   well_tmp = meta_tmp$well[i]
  #   chain_tmp = meta_tmp$chain[i]
  #   #print(file_full_tmp)
  #   if(file_type == ".tsv") {
  #     df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
  #       mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
  #       mutate(well_chain = paste(well, chain))
  #   } else if(file_type == ".parquet") {
  #     df_tmp = arrow::read_parquet(file_full_tmp) %>%
  #       mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
  #       mutate(well_chain = paste(well, chain))
  #     if(!is.null(columns)) df_tmp = df_tmp[,columns, drop = FALSE]
  #   }
  #   return(df_tmp)
  # }, mc.cores = as.integer(nproc), mc.preschedule = FALSE) %>% bind_rows()  %>%
  #   dplyr::rename(j = allJHitsWithScore, v = allVHitsWithScore) %>%
  #   mutate(j = gsub("\\(.*", "", j), v = gsub("\\(.*", "", v))

  #ll = well_data_df_to_sparse(df)
  return(df)

}


well_data_df_to_sparse = function(df, cols, transpose = F) {
  #all_wells = get_well_subset()
  #cols = all_wells[all_wells %in% unique(df$well)]
  col_df = tibble(well = cols, col_id = 1:length(cols))
  df$row_id <- vctrs::vec_group_id(df[, c("targetSequences", "v", "j")])
  df_row_meta = df[!duplicated(row_id),c("targetSequences", "v", "j", "aaSeqCDR3", "chain", "row_id")]
  df$col_id = col_df$col_id[match(df$well, col_df$well)]
  if(!transpose) {
    mat = sparseMatrix(
      i = df$row_id,
      j = df$col_id,
      x = df$readCount,
      dims = c(max(df$row_id), max(df$col_id))
    )
  } else {
    mat = sparseMatrix(
      j = df$row_id,
      i = df$col_id,
      x = df$readCount,
      dims = c(max(df$col_id), max(df$row_id))
    )
  }
  ll = list(row_data = df_row_meta, col_data = col_df, matrix = mat)
  return(ll)
}


well_data_df_to_sparse_multi = function(df_list) {
  all_wells = get_well_subset()
  cols_list = lapply(df_list, function(df) all_wells[all_wells %in% unique(df$well)])
  col_df_list = lapply(cols_list, function(cols) tibble(well = cols, col_id = 1:length(cols)))
  df = bind_rows(df_list) %>%
    dplyr::rename(j = allJHitsWithScore, v = allVHitsWithScore) %>%
    mutate(j = gsub("\\(.*", "", j), v = gsub("\\(.*", "", v))
  #cols = all_wells[all_wells %in% unique(df$well)]
  #col_df = tibble(well = cols, col_id = 1:length(cols))
  df$row_id <- vctrs::vec_group_id(df[, c("targetSequences", "v", "j")])
  df_row_meta = df[!duplicated(row_id),c("targetSequences", "v", "j", "aaSeqCDR3", "chain", "row_id")]
  nrows = nrow(df_row_meta)

  df_list = split(df, df$sample_num)
  mat_list = lapply(df_list, function(df_tmp) {
    i = unique(df_tmp$sample_num)
    col_df = col_df_list[[i]]
    df_tmp$col_id = col_df$col_id[match(df_tmp$well, col_df$well)]
    mat = Matrix::sparseMatrix(
      i = df_tmp$row_id,
      j = df_tmp$col_id,
      x = df_tmp$readCount,
      dims = c(nrows, max(df_tmp$col_id))
    )
    return(mat)
  })
  ll = list(row_data = df_row_meta, col_data_list = col_df_list, matrix_list = mat_list)
  return(ll)
}


save_well_counts_h5 = function(matrix, prefix, group = "matrix") {
  HDF5Array::writeTENxMatrix(matrix, filepath = paste(prefix,".h5",sep=""), group = "matrix")
  return(invisible(NULL))
}

load_wells_counts_h5 = function(prefix, folder) {
  alpha_mat_file = file.path(folder, paste(prefix, "alpha.h5", sep="_"))
  beta_mat_file = file.path(folder, paste(prefix, "beta.h5", sep="_"))

  alpha_row_file = file.path(folder, paste(prefix, "alpha.parquet", sep="_"))
  beta_row_file = file.path(folder, paste(prefix, "beta.parquet", sep="_"))

  col_file = file.path(folder, paste(prefix, "cols.parquet", sep="_"))

  alpha_mat = HDF5Array::TENxMatrix(alpha_mat_file, group = "matrix") %>% t() %>% as("dgCMatrix")
  beta_mat = HDF5Array::TENxMatrix(beta_mat_file, group = "matrix") %>% t() %>% as("dgCMatrix")

  alpha_row = arrow::read_parquet(alpha_row_file)
  beta_row = arrow::read_parquet(beta_row_file)
  cols = arrow::read_parquet(col_file)

  out_list = list(alpha = alpha_mat, beta = beta_mat, rows_alpha = alpha_row, rows_beta = beta_row, cols = cols)
  return(out_list)
}

lazy_load_parquet = function(file) {
  con <- duckdb::dbConnect(duckdb::duckdb())
  call = paste("read_parquet(", "'", file, "'", ")", sep = "")
  #tbl <- tbl(con, "read_parquet('alpha_test.parquet')")
  df <- dplyr::tbl(con, call)
  return(df)
}

load_wells_counts_rds = function(folder, prefix, lazy = TRUE) {
  alpha_mat_file = file.path(folder, paste(prefix, "alpha.rds", sep="_"))
  beta_mat_file = file.path(folder, paste(prefix, "beta.rds", sep="_"))

  alpha_row_file = file.path(folder, paste(prefix, "alpha.parquet", sep="_"))
  beta_row_file = file.path(folder, paste(prefix, "beta.parquet", sep="_"))

  col_file = file.path(folder, paste(prefix, "cols.parquet", sep="_"))

  alpha_mat = readRDS(alpha_mat_file)
  beta_mat = readRDS(beta_mat_file)

  if(lazy) {
    alpha_row = lazy_load_parquet(alpha_row_file)
    beta_row = lazy_load_parquet(beta_row_file)
    cols = arrow::read_parquet(col_file)
  } else {
    alpha_row = arrow::read_parquet(alpha_row_file)
    beta_row = arrow::read_parquet(beta_row_file)
    cols = arrow::read_parquet(col_file)
  }

  out_list = list(alpha = alpha_mat, beta = beta_mat,
                  alpha_meta = alpha_row, beta_meta = beta_row, col_meta = cols,
                  alpha_meta_file = alpha_row_file, beta_meta_file = beta_row_file
                  )
  return(out_list)
}

write_well_data_to_binary = function(
    folder_in,
    folder_out,
    prefix,
    wells = get_well_subset(1:16,1:24),
    in_file_type = c(".tsv", ".parquet"),
    out_file_type = c(".rds", ".h5"),
    well_pos=3,
    parallel = FALSE,
    nproc = data.table::getDTthreads(),
    columns = NULL,
    max_files = Inf,
    periods_to_underscores = TRUE,
    to_sparse_matrix = TRUE
    ) {

  .create_folder(folder_out)
  df_alpha = load_well_data_to_df_single(
    folder_in,
    wells = wells,
    chain = "alpha",
    file_type = in_file_type,
    well_pos = well_pos,
    parallel = parallel,
    nproc = nproc,
    columns = columns,
    max_files = max_files,
    periods_to_underscores = periods_to_underscores
  )
  df_beta = load_well_data_to_df_single(
    folder_in,
    wells = wells,
    chain = "beta",
    file_type = in_file_type,
    well_pos = well_pos,
    parallel = parallel,
    nproc = nproc,
    columns = columns,
    max_files = max_files,
    periods_to_underscores = periods_to_underscores
  )
  if(to_sparse_matrix) {
    ll_alpha = well_data_df_to_sparse(df_alpha, cols = wells)
    ll_beta = well_data_df_to_sparse(df_beta, cols = wells)
    check = all.equal(ll_alpha$col_data, ll_beta$col_data)
    if(!check) stop("well ordering for alpha and beta doesn't match")
    if(!out_file_type %in% c(".rds", ".h5")) stop("'out_file_type' must be either '.rds' or '.h5'")
    if(out_file_type == ".h5") {
      HDF5Array::writeTENxMatrix(
        ll_alpha$matrix,
        filepath = file.path(folder_out, paste(prefix,"alpha.h5",sep="_")),
        group = "matrix")
      HDF5Array::writeTENxMatrix(
        ll_beta$matrix,
        filepath = file.path(folder_out, paste(prefix,"beta.h5",sep="_")),
        group = "matrix")
    } else if(out_file_type == ".rds") {
      saveRDS(ll_alpha$matrix, file = file.path(folder_out, paste(prefix,"alpha.rds",sep="_")))
      saveRDS(ll_beta$matrix, file = file.path(folder_out, paste(prefix,"beta.rds",sep="_")))
    }
    arrow::write_parquet(
      ll_alpha$row_data,
      file.path(folder_out, paste(prefix,"alpha.parquet",sep="_")))
    arrow::write_parquet(
      ll_beta$row_data,
      file.path(folder_out, paste(prefix,"beta.parquet",sep="_")))
    arrow::write_parquet(
      ll_beta$col_data,
      file.path(folder_out, paste(prefix,"cols.parquet",sep="_")))
  } else {
    arrow::write_parquet(
      df_alpha,
      file.path(folder_out, paste(prefix,"alpha.parquet",sep="_")))
    arrow::write_parquet(
      df_beta,
      file.path(folder_out, paste(prefix,"beta.parquet",sep="_")))
  }
  return(invisible(NULL))
}
