#' Install python environment
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function install a python environment via the `basilisk` package that is
#' needed for some functions: \code{\link{TCRdist}()}, \code{\link{cluster_tcrs}()}, and \code{\link{run_pairing}()}.
#'
#' It will run automatically when needed, but is available to the user if needed.
#'
#' @family python-deps
#' @export
#'

install_python_env <- function() {
  message("Installing/updating the Python environment for TIRTLtools...")
  install_python_env_silently()
  message("Done.")
  return(invisible(NULL))
}
