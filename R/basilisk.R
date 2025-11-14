.TIRTLtools_deps <- c(
  "python==3.12.9",
  "numpy==2.3.0",
  "scipy==1.15.2",
  "pandas==2.2.2",
  "mlx==0.23.0"#,
  #"cupy==13.4.1"
)

TIRTLtools_env <- basilisk::BasiliskEnvironment(
  envname = "TIRTLtools_env",
  pkgname = "TIRTLtools",
  packages = .TIRTLtools_deps
)
