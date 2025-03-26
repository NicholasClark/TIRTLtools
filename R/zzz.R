.onLoad <- function(...) {
  #reticulate::py_require(c("mlx", "cupy", "numpy", "pandas", "datetime", "sys"))
  reticulate::use_condaenv("nick_main")
  data("params", "submat", package="TIRTLtools", envir=parent.env(environment()))
  #np <- import("numpy", delay_load = TRUE)
  utils = reticulate::import_from_path("utils", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
  TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
}
