# .env_installed_flag <- function() {
#   # Per-user cache dir for TIRTLtools
#   file.path(tools::R_user_dir("TIRTLtools", which = "cache"),
#             "env_installed.flag")
# }



# mark_env_installed <- function() {
#   dir.create(dirname(.env_installed_flag()),
#              showWarnings = FALSE, recursive = TRUE)
#   file.create(.env_installed_flag())
# }

# Internal: trigger basilisk env creation without doing real work
install_python_env_silently <- function() {
  basilisk::basiliskRun(
    env = TIRTLtools_py_env,
    fun = function() NULL
  )
  #mark_env_installed()
  invisible(TRUE)
}


ensure_python_env <- function(ask = TRUE) {

  if (env_is_installed()) {
    return(invisible(TRUE))
  }

  # Non-interactive or explicit "don't ask": auto-install
  if (!interactive() || !isTRUE(ask)) {
    install_python_env_silently()
    return(invisible(TRUE))
  }

  # Interactive + ask = TRUE: prompt the user
  version = packageVersion("TIRTLtools")
  message("This function requires installation of a Python environment.")
  message(paste("This is the first time you are using the Python backend for TIRTLtools v", version, sep = ""))


  ans <- utils::menu(
    choices = c("Yes, install now", "No, cancel"),
    title = "Allow TIRTLtools to install its Python environment now?"
  )

  if (ans != 1) {
    stop("Installation of the Python environment was cancelled by the user.")
  }

  install_python_env_silently()
  invisible(TRUE)
}
