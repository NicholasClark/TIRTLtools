base_deps <- c(
  "python==3.12.9",
  "numpy==2.3.0",
  "scipy==1.15.2",
  "pandas==2.2.2"
)

silicon_deps <- c(base_deps,"mlx==0.23.0")
nvidia_deps <- c(base_deps,"cupy==13.4.1")

TIRTLtools_env_no_gpu <- basilisk::BasiliskEnvironment(
  envname = "TIRTLtools_env_no_gpu",
  pkgname = "TIRTLtools",
  packages = base_deps
)

TIRTLtools_env_nvidia <- basilisk::BasiliskEnvironment(
  envname = "TIRTLtools_env",
  pkgname = "TIRTLtools",
  packages = nvidia_deps
)

TIRTLtools_env_silicon <- basilisk::BasiliskEnvironment(
  envname = "TIRTLtools_env",
  pkgname = "TIRTLtools",
  packages = silicon_deps
)
