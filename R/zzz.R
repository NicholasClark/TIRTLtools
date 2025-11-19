# .onAttach <- function(libname, pkgname) {
#   if (interactive()) {
#     packageStartupMessage(
#       "TIRTLtools: to pre-install the Python environment, run:\n",
#       "  TIRTLtools::install_python_env()"
#     )
#   }
# }




# .onLoad <- function(...) {
#   ### code to install python dependencies (conda) ----------
#   install_python_deps_conda <- function(envname = "TIRTLtools") {
#     # Check if Python is available
#     # if (!reticulate::py_available(initialize = FALSE)) {
#     #   stop("Python is not available. Please install Python before continuing.")
#     # }
#     if (!reticulate::condaenv_exists(envname)) {
#       reticulate::conda_create(envname)
#     }
#     reticulate::use_condaenv(envname)
#
#     # non-gpu packages needed
#     #py_packages = c("numpy", "scipy", "pandas", "time", "os", "subprocess", "platform", "importlib")
#     py_packages = c("numpy", "scipy", "pandas")
#     install_pkgs <- py_packages[!vapply(py_packages, reticulate::py_module_available, logical(1))]
#
#     # Detect OS and GPU support
#     sys_info <- Sys.info()
#     is_macos <- sys_info["sysname"] == "Darwin"
#     #is_apple_silicon <- is_macos && grepl("arm64", R.version$arch)
#     is_apple_silicon <- is_macos && grepl("arm64", sys_info['machine'])
#
#
#     # macOS Apple Silicon (M1/M2/M3 chips) → try mlx
#     if (is_apple_silicon) {
#       if (! reticulate::py_module_available("mlx") ) {
#         message("Apple Silicon detected. Will attempt to install 'mlx'.")
#         install_pkgs = c(install_pkgs, "mlx")
#       }
#     }
#
#     # Try to detect Nvidia GPU
#     has_nvidia_gpu <- FALSE
#     if (!is_macos) {
#       # Try system call to 'nvidia-smi'
#       nvidia_smi <- try(suppressWarnings(system("nvidia-smi", intern = TRUE, ignore.stderr = TRUE)))
#       if (length(nvidia_smi) > 0 && class(nvidia_smi) != "try-error") {
#         has_nvidia_gpu <- TRUE
#         if (! reticulate::py_module_available("cupy") ) {
#           message("NVIDIA GPU detected. Will attempt to install 'cupy'.")
#           install_pkgs = c(install_pkgs, "cupy")
#         }
#       }
#     }
#
#     # Install the selected Python package
#     if( length(install_pkgs) > 0 ) {
#       reticulate::py_install(
#         packages = install_pkgs,
#         envname = envname,
#         method = "auto",
#         #method = "conda",
#         pip = TRUE
#       )
#     }
#   }
#
#   install_python_deps_no_conda <- function(envname = "TIRTLtools") {
#     reticulate::py_available(initialize = TRUE)
#     # non-gpu packages needed
#     #py_packages = c("numpy", "scipy", "pandas", "time", "os", "subprocess", "platform", "importlib")
#     py_packages = c("numpy", "scipy", "pandas")
#     install_pkgs <- py_packages[!vapply(py_packages, reticulate::py_module_available, logical(1))]
#
#     # Detect OS and GPU support
#     sys_info <- Sys.info()
#     is_macos <- sys_info["sysname"] == "Darwin"
#     #is_apple_silicon <- is_macos && grepl("arm64", R.version$arch)
#     is_apple_silicon <- is_macos && grepl("arm64", sys_info['machine'])
#
#
#     # macOS Apple Silicon (M1/M2/M3 chips) → try mlx
#     if (is_apple_silicon) {
#       if (! reticulate::py_module_available("mlx") ) {
#         #message("Apple Silicon detected. Will attempt to install 'mlx'.")
#         install_pkgs = c(install_pkgs, "mlx")
#       }
#     }
#
#     # Try to detect Nvidia GPU
#     has_nvidia_gpu <- FALSE
#     if (!is_macos) {
#       # Try system call to 'nvidia-smi'
#       nvidia_smi <- try(suppressWarnings(system("nvidia-smi", intern = TRUE, ignore.stderr = TRUE)))
#       if (length(nvidia_smi) > 0 && class(nvidia_smi) != "try-error") {
#         has_nvidia_gpu <- TRUE
#         if (! reticulate::py_module_available("cupy") ) {
#           #message("NVIDIA GPU detected. Will attempt to install 'cupy'.")
#           install_pkgs = c(install_pkgs, "cupy")
#         }
#       }
#     }
#
#     # Install the selected Python package
#     if( length(install_pkgs) > 0 ) {
#     #   reticulate::py_install(
#     #     packages = install_pkgs,
#     #     #envname = envname,
#     #     method = "auto",
#     #     #method = "conda",
#     #     pip = TRUE
#     #   )
#       reticulate::py_require(install_pkgs)
#     }
#   }
#
#   run_conda = function() {
#     install_python_deps_conda()
#     reticulate::use_condaenv("TIRTLtools")
#     reticulate::py_available(initialize = TRUE)
#   }
#
#   # install_python_deps_alt = function() {
#   #   reticulate::py_available(initialize = TRUE)
#   #   reticulate::py_require(c("numpy", "scipy", "pandas", "time", "os", "subprocess", "platform", "importlib"))
#   # }
#
#   #run_conda()
#   #install_python_deps_alt()
#   install_python_deps_no_conda()
#
#   data("params", "submat", package="TIRTLtools", envir=parent.env(environment()))
#
#   #### previous code (just loads conda environment) --------
#   #reticulate::py_require(c("mlx", "cupy", "numpy", "pandas", "datetime", "sys"))
#   #reticulate::use_condaenv("nick_main")
#   # reticulate::use_condaenv("TIRTLtools")
#   # reticulate::py_available(initialize = TRUE)
#   # data("params", "submat", package="TIRTLtools", envir=parent.env(environment()))
#   #np <- import("numpy", delay_load = TRUE)
#   # utils = reticulate::import_from_path("utils", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
#   # TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
#   #utils = reticulate::import_from_path("utils", path = "inst/python", convert = TRUE, delay_load = TRUE)
#   #TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = "inst/python", convert = TRUE, delay_load = TRUE)
#   ######
#   return(NULL)
# }
