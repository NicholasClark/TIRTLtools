### function for individual well files from a folder

# test = read_well_data(folder_path = "/Users/nclark52/git/TIRTL_nick/data/exp3_clones/TCR_clones_ID03/", wells = get_well_subset(1:2,1:2), chain = "beta")

read_well_data = function(folder_path, wells = get_well_subset(1:16,1:24),
                          well_pos=3, chain = c("both", "alpha", "beta"), nproc = 1L,
                          columns = NULL,
                          max_files = Inf) {
  chain = chain[1]

  files = list.files(path = folder_path, full.names = F)
  files_full = file.path(folder_path, files)
  files_no_ending = gsub("\\.tsv", "", files)
  #files_no_ending = gsub("\\.", "_", files_no_ending) ## replace periods with underscores if necessary

  wells_all = sapply(strsplit(files_no_ending, "_"), function(x) x[[well_pos]])
  ## make sure wells in the string are 2-3 characters long and only alphanumeric -- necessary for some data
  wells_all = sapply(wells_all, function(x) substr(x, 1, 3))
  wells_all = gsub("[^[:alnum:] ]", "", wells_all)

  chains_all = sapply(strsplit(files_no_ending, "_"), function(x) rev(x)[[1]])
  chains_all[chains_all=="TRA"] = "alpha"
  chains_all[chains_all=="TRB"] = "beta"

  wells_keep = wells_all %in% wells
  if(chain %in% c("alpha","beta")) {
    chains_keep = chains_all == chain
  } else {
    chains_keep = TRUE
  }
  keep = which(wells_keep & chains_keep)

  meta_tmp = tibble(well = wells_all[keep], chain = chains_all[keep], file_short = files[keep], file_full = files_full[keep])
  print(paste(dim(meta_tmp)[1], "files total"))
  n_files = min(max_files, dim(meta_tmp)[1])
  print(paste("Loading", n_files, "files"))
  ll = parallel::mclapply(1:n_files, function(i) {
    print(i)
    file_full_tmp = meta_tmp$file_full[i]
    file_short_tmp = meta_tmp$file_short[i]
    well_tmp = meta_tmp$well[i]
    chain_tmp = meta_tmp$chain[i]
    #print(file_full_tmp)
    df_tmp = data.table::fread(file_full_tmp, select = columns) %>%
      mutate(chain = chain_tmp, well = well_tmp, file_short = file_short_tmp) %>%
      mutate(well_chain = paste(well, chain))
    return(df_tmp)
    }, mc.cores = as.integer(nproc))
  #names(ll) = files_no_ending[keep]
  return(ll)
}
