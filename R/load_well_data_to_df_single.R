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
      df_tmp = nanoparquet::read_parquet(file_full_tmp) %>%
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
  #     df_tmp = nanoparquet::read_parquet(file_full_tmp) %>%
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
