#' Install python environment
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function install a python environment via the `basilisk` package that is
#' needed for some functions: \code{\link{TCRdist}()}, \code{\link{cluster_tcrs}()}, and \code{\link{run_pairing}()}.
#'
#' It will run automatically when needed, but is available to the user if needed.
#'
#' @param force whether to force installation (default if FALSE)
#'
#' @return returns NULL
#'
#' @family python-deps
#'

install_python_env <- function(force = FALSE) {
  message("Installing/updating the Python environment for TIRTLtools...")
  install_python_env_silently(force = force)
  message("Done.")
  return(invisible(NULL))
}
