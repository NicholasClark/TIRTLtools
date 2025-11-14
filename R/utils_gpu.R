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

.choose_basilisk_env <- function(backend = c("auto", "cpu", "cupy", "mlx")) {
  backend <- match.arg(backend)

  if (backend == "cpu")  return(TIRTLtools_env_no_gpu)
  if (backend == "cupy") return(TIRTLtools_env_nvidia)
  if (backend == "mlx")  return(TIRTLtools_env_silicon)

  ## backend == "auto"
  if (.is_apple_silicon()) {
    return(TIRTLtools_env_silicon)
  }

  if (.has_nvidia_gpu()) {
    return(TIRTLtools_env_nvidia)
  }

  return(TIRTLtools_env_no_gpu)
}
