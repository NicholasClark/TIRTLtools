download_dataset = function(dataset = "SJTRC_TIRTLseq_minimal") {

  file_info = .get_dataset_info(dataset)
  dir.create(dataset, showWarnings = FALSE)

  for (i in 1:dim(file_info)[1]) {
    url = file_info$download_url[i]
    file_name = file_info$name[i]
    msg = paste("Downloading file:", file_name, " to directory: ", dataset) %>% .add_newline()
    cat(msg)
    dest = file.path(dataset, basename(url))
    download.file(url, destfile = dest, mode = "wb")
  }
  return(invisible(NULL))
}

.get_dataset_info = function(dataset = "SJTRC_TIRTLseq_minimal") {
  ds_options = c("SJTRC_TIRTLseq_minimal")
  dataset = dataset[1]
  if(!dataset %in% ds_options) stop("dataset not available")

  minimal_data_df = data.frame(name = c("cd8_tp1_v2_TIRTLoutput.tsv.gz", "cd8_tp1_v2_pseudobulk_TRA.tsv.gz",
                                        "cd8_tp1_v2_pseudobulk_TRB.tsv.gz", "cd8_tp2_v2_TIRTLoutput.tsv.gz",
                                        "cd8_tp2_v2_pseudobulk_TRA.tsv.gz", "cd8_tp2_v2_pseudobulk_TRB.tsv.gz"
  ), download_url = c("https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp1_v2_TIRTLoutput.tsv.gz",
                      "https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp1_v2_pseudobulk_TRA.tsv.gz",
                      "https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp1_v2_pseudobulk_TRB.tsv.gz",
                      "https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp2_v2_TIRTLoutput.tsv.gz",
                      "https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp2_v2_pseudobulk_TRA.tsv.gz",
                      "https://raw.githubusercontent.com/NicholasClark/TIRTLtools_data/main/SJTRC_minimal_dataset/cd8_tp2_v2_pseudobulk_TRB.tsv.gz"
  ))

  if(dataset == "SJTRC_TIRTLseq_minimal") file_info = minimal_data_df
  return(file_info)
}

# .get_dataset_info_gh = function(dataset = "minimal") {
#   ds_options = c("minimal")
#   dataset = dataset[1]
#   if(!dataset %in% ds_options) stop("dataset not available")
#   if(dataset == "minimal") {
#     files = gh::gh("GET /repos/:owner/:repo/contents/:path",
#                owner = "NicholasClark",
#                repo = "TIRTLtools_data",
#                path = "SJTRC_minimal_dataset")
#
#     file_info = lapply(files, function(x) list(name = x$name, download_url = x$download_url))
#     file_info = do.call(rbind, lapply(file_info, as.data.frame))
#     return(file_info)
#   }
# }
