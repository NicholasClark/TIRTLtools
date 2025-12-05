## original deps -- working for basilisk 1.20
# base_deps <- c(
#   "python==3.12.9",
#   "numpy==2.3.0",
#   "scipy==1.15.2",
#   "pandas==2.2.2"
# )

## deps for basilisk 1.22 (using pyenv)
base_deps <- c(
  #"python==3.14",
  "numpy==2.3.5",
  "scipy==1.16.3",
  "pandas==2.3.3"
)

#silicon_deps <- c(base_deps,"mlx==0.23.0") ## working for basilisk 1.20
#nvidia_deps <- c(base_deps,"cupy==13.4.1") ## working for basilisk 1.20

.has_nvidia_gpu <- function() {
  # Simple heuristic: uses nvidia-smi or driver file
  if (.Platform$OS.type != "unix") return(FALSE)

  if (nzchar(Sys.which("nvidia-smi"))) {
    # Check that it actually works
    status <- suppressWarnings(
      system("nvidia-smi -L > /dev/null 2>&1")
    )
    if (identical(status, 0L)) return(TRUE)
  }

  # Linux driver file
  if (file.exists("/proc/driver/nvidia/version")) return(TRUE)

  FALSE
}

.is_apple_silicon <- function() {
  si <- Sys.info()
  identical(si[["sysname"]], "Darwin") && grepl("arm64", si[["machine"]])
}

if(.has_nvidia_gpu()) {
  base_deps = c(base_deps,"cupy==13.6.0")
}
if (.is_apple_silicon()) {
  base_deps <- c(base_deps,"mlx==0.29.4")
}

envname = "TIRTLtools_py_env"
pkgname = "TIRTLtools"

TIRTLtools_py_env <- basilisk::BasiliskEnvironment(
  envname = envname,
  pkgname = pkgname,
  packages = base_deps
)

get_py_deps = function() {
  return(base_deps)
}

env_is_installed <- function() {
  if("getExternalDir" %in% getNamespaceExports("basilisk")) {
    exdir <- basilisk::getExternalDir()
  } else if("getExternalDir" %in% getNamespaceExports("basilisk.utils")) {
    exdir <- basilisk.utils::getExternalDir()
  } else {
    return(TRUE) ## not sure if python installed or not, don't ask to install
  }
  envdir <- file.path(exdir, pkgname, packageVersion(pkgname))
  envpath <- file.path(envdir, envname) ## envname defined above this function
  return(file.exists(envpath))
  # if (is.null(pkgname)) {
  #   envpath <- envname
  # } else if (basilisk::useSystemDir()) {
  #   envpath <- file.path(basilisk:::.get_env_system_dir(pkgname, installed = TRUE),
  #                        envname)
  # } else {
  #   exdir <- basilisk::getExternalDir()
  #   envdir <- file.path(exdir, pkgname, packageVersion(pkgname))
  #   envpath <- file.path(envdir, envname)
  # }
  # file.exists(envpath)
  #file.exists(.env_installed_flag())
}

# TIRTLtools_env_no_gpu <- basilisk::BasiliskEnvironment(
#   envname = "TIRTLtools_env_no_gpu",
#   pkgname = "TIRTLtools",
#   packages = base_deps
# )
#
# TIRTLtools_env_nvidia <- basilisk::BasiliskEnvironment(
#   envname = "TIRTLtools_env",
#   pkgname = "TIRTLtools",
#   packages = nvidia_deps
# )
#
# TIRTLtools_env_silicon <- basilisk::BasiliskEnvironment(
#   envname = "TIRTLtools_env",
#   pkgname = "TIRTLtools",
#   packages = silicon_deps
# )

# deps <- c(
#   "python==3.12.9",
#   "numpy==2.3.0",
#   "scipy==1.15.2",
#   "pandas==2.2.2",
#   "mlx==0.23.0",
#   "cupy==13.4.1"
# )
#
# TIRTLtools_env_all <- basilisk::BasiliskEnvironment(
#   envname = "TIRTLtools_env_all",
#   pkgname = "TIRTLtools",
#   packages = deps
# )
