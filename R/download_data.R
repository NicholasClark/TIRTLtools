download_data = function(dataset = c("SJTRC_minimal.qs2", "SJTRC_longitudinal.qs2", "exp3_tp1_cd8.tar.gz")) {
  dataset = dataset[1]
  cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
  assert_choice(dataset, choices = c("SJTRC_minimal.qs2", "SJTRC_longitudinal.qs2", "exp3_tp1_cd8.tar.gz"))
  
  tag = "data-v1"
  repo = "NicholasClark/TIRTLtools"
  asset = dataset
  size_df = piggyback:::pb_list(repo = repo, tag = tag) %>% 
    filter(file_name == dataset)
  file_size = scales::label_bytes()(ifelse(nrow(size_df) >= 1, size_df$size[1], NA))
  download_confirm = .confirm(glue("Download {dataset} ({file_size})?"))
  if(!download_confirm) {
    message("Download cancelled.")
    return(invisible(FALSE))
  }
  tag_dir <- file.path(cache_root, tag)
  dir.create(tag_dir, recursive = TRUE, showWarnings = FALSE)
  file_path <- file.path(tag_dir, asset)
  piggyback::pb_download(
    file = asset,
    repo = repo,
    tag = tag,
    dest = tag_dir,
    overwrite = TRUE
  )
  msg = paste("Dataset downloaded to:", "\n", file_path)
  message(msg)
  if(grepl("\\.tar\\.gz$", file_path)) extract_data(file_path)
  return(invisible(TRUE))
}

## note: to work correctly, "{dataset}.tar.gz" needs to extract to "parent_dir/{dataset}/" directory
extract_data = function(path) {
  parent_dir = dirname(path)
  extract_dir = file.path(parent_dir, sub("\\.tar\\.gz", "", basename(path)))
  message("Extracting archive...")
  utils::untar(path, exdir = parent_dir)
  files = dir(extract_dir)
  n_files = length(files)
  msg = paste("Extracted", n_files, "files to:", "\n", extract_dir)
  message(msg)
  msg = paste(c( head(files,3), "..."), collapse = "\n")
  message(msg)
  return(invisible(NULL))
}

clean_cache = function() {
  unlink(file.path(tools::R_user_dir("TIRTLtools", which = "cache"), "*"), recursive = TRUE)
  return(invisible(NULL))
}

# download_well_data = function(dataset = "exp3_tp1_cd8") {
#   dataset = dataset[1]
#   cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
#   assert_choice(dataset, choices = "exp3_tp1_cd8")

#   if(dataset == "exp3_tp1_cd8") {
#     tag = "data-v1"
#     repo = "NicholasClark/TIRTLtools"
#     asset = "exp3_tp1_cd8.tar.gz"

#     tag_dir <- file.path(cache_root, tag)
#     dir.create(tag_dir, recursive = TRUE, showWarnings = FALSE)
#     final_path = file.path(tag_dir, "exp3_tp1_cd8")
#     check = length(dir(final_path)) == 384
#     # if(check && !force) {
#     #   msg = paste("Well-level data files already downloaded to:", "\n", final_path)
#     #   message(msg)
#     #   return(final_path)
#     # }


#     archive_path <- file.path(tag_dir, asset)
#     #if(file.exists(archive_path)) unlink(archive_path) ## cleaning up if needed
#     #if(file.exists(final_path)) unlink(final_path, recursive = TRUE) ## cleaning up if needed
#     size_df = piggyback:::pb_list(repo = repo, tag = tag) %>% 
#       filter(file_name == asset)
#     file_size = scales::label_bytes()(ifelse(nrow(size_df) >= 1, size_df$size[1], NA))
#     .confirm(glue("Download {asset} ({file_size})?"))
#     piggyback::pb_download(
#       file = asset,
#       repo = repo,
#       tag = tag,
#       dest = tag_dir,
#       overwrite = TRUE
#     )
#     #download_url = glue::glue("https://github.com/{repo}/releases/download/{tag}/{asset}")
#     #curl::multi_download(download_url, file.path(tag_dir, asset), resume = TRUE, progress = TRUE)
#     #curl::curl_download(download_url, file.path(tag_dir, asset), quiet = FALSE)

#     message("Extracting archive...")
#     utils::untar(archive_path, exdir = tag_dir)
#     files = dir(final_path)
#     n_files = length(files)
#     msg = paste("Extracted", n_files, "files:")
#     message(msg)
#     msg = paste(c( head(files,3), "..."), collapse = "\n")
#     message(msg)
#     #file.remove(archive_path)

#     msg = paste("Well-level data files downloaded to:", "\n", final_path)
#     message(msg)
#     return(final_path)
#   }
# }