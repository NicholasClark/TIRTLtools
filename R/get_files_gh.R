get_files_gh = function(user, repo, folder) {
  api_req = glue::glue("GET /repos/{user}/{repo}/contents/{folder}")
  files = sapply(gh::gh(api_req), function(x) x$name)
  url_base = glue::glue("https://github.com/{user}/{repo}/raw/refs/heads/main/{folder}/")
  urls = paste(url_base, files, sep = "")
  return(urls)
}


get_files_gh_SJTRC_minimal = function() {
  get_files_gh("NicholasClark", "TIRTLtools_data", "SJTRC_minimal_dataset")
}

