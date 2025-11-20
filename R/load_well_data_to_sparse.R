
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
    nproc = 1L,
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

load_well_data_to_sparse_single = function(folder_path, wells = get_well_subset(1:16,1:24),
                                       well_pos=3, chain = c("alpha", "beta"), nproc = 1L,
                                       columns = NULL,
                                       max_files = Inf, periods_to_underscores = TRUE) {
  args <- as.list(environment())
  df = do.call("load_well_data_to_df_single", args)
  ll = well_data_df_to_sparse(df)
  return(ll)
}

load_well_data_to_df_single = function(folder_path, wells = get_well_subset(1:16,1:24),
                                    well_pos=3, chain = c("alpha", "beta"), nproc = 1L,
                                    columns = NULL,
                                    max_files = Inf, periods_to_underscores = TRUE) {
  chain = chain[1]
  if(!chain %in% c("alpha", "beta")) stop("'chain' must be 'alpha' or 'beta'.")
  chk1 = list.files(path = folder_path, pattern = "\\.gz$")
  if(length(chk1) > 0) warning("Folder contains .gz files -- check to make sure they have all been extracted to .tsv")
  if(chain == "alpha") files = list.files(path = folder_path, full.names = F, pattern = "_TRA\\.tsv$")
  if(chain == "beta") files = list.files(path = folder_path, full.names = F, pattern = "_TRB\\.tsv$")
  files_full = file.path(folder_path, files)
  files_no_ending = gsub("\\.tsv", "", files)
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
  df = parallel::mclapply(1:n_files, function(i) {
    cat("|")
    file_full_tmp = meta_tmp$file_full[i]
    file_short_tmp = meta_tmp$file_short[i]
    well_tmp = meta_tmp$well[i]
    chain_tmp = meta_tmp$chain[i]
    #print(file_full_tmp)
    df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
      mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
      mutate(well_chain = paste(well, chain))
    return(df_tmp)
  }, mc.cores = as.integer(nproc)) %>% bind_rows()  %>%
    dplyr::rename(j = allJHitsWithScore, v = allVHitsWithScore) %>%
    mutate(j = gsub("\\(.*", "", j), v = gsub("\\(.*", "", v))
  #ll = well_data_df_to_sparse(df)
  return(df)
}


well_data_df_to_sparse = function(df) {
  all_wells = get_well_subset()
  cols = all_wells[all_wells %in% unique(df$well)]
  col_df = tibble(well = cols, col_id = 1:length(cols))
  df$row_id <- vctrs::vec_group_id(df[, c("targetSequences", "v", "j")])
  df_row_meta = df[!duplicated(row_id),c("targetSequences", "v", "j", "aaSeqCDR3", "chain", "row_id")]
  df$col_id = col_df$col_id[match(df$well, col_df$well)]
  mat = sparseMatrix(
    i = df$row_id,
    j = df$col_id,
    x = df$readCount,
    dims = c(max(df$row_id), max(df$col_id))
  )
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
