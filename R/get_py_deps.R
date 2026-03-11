get_py_deps = function(with_versions = FALSE) {
  base_deps <- c(
    #"python==3.14",
    "numpy==2.3.5",
    "scipy==1.16.3",
    "pandas==2.3.3"
  )
  if(.has_nvidia_gpu()) {
    base_deps = c(base_deps,"cupy==13.6.0")
  }
  if (.is_apple_silicon()) {
    base_deps <- c(base_deps,"mlx==0.29.4")
  }

  if(with_versions) return(base_deps)
  return( gsub("==.*", "", base_deps) )
}

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
