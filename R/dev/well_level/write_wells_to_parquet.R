write_wells_to_parquet = function(folder_in, folder_out, nproc = data.table::getDTthreads(), max_files = Inf) {
  tictoc::tic()
  .create_folder(folder_out)
  files_tsv = dir(folder_in, pattern = "\\.tsv$", full.names = TRUE)
  n_files = min(c(length(files_tsv), max_files))

  parallel::mclapply(1:n_files, function(i) {
    cat(paste(i, ""))
    file = files_tsv[i]
    tmp = data.table::fread(file)
    file_out_base = gsub("\\.tsv$",".parquet", basename(file))
    nanoparquet::write_parquet(tmp, file.path(folder_out, file_out_base))
    return(NULL)
  }, mc.cores = nproc)

  msg = paste(n_files, ".parquet files written to folder:", folder_out)
  message(msg)
  tictoc::toc()
  return(invisible(NULL))
}
