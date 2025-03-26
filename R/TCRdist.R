## TCRdist main function

# tcr_all = arrow::read_parquet("../TCRdist_gpu/TCRdist_gpu/data/vdj_and_covid_032425.parquet")
# tcr1 = tcr_all[1:20000,]
# tcr2 = tcr_all[1:20000,]
# tmp = TCRdist(tcr1 = tcr1, tcr2 = tcr2)
TCRdist = function(tcr1=NULL, tcr2=NULL, params_df = NULL, submat = NULL, tcrdist_cutoff=90, chunk_size=1000, print_chunk_size=1000, print_res = TRUE) {
  reticulate::use_condaenv("nick_main")
  #np <- import("numpy", delay_load = TRUE)
  utils = reticulate::import_from_path("utils", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
  TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = "../TCRdist_gpu/TCRdist_gpu/", convert = TRUE, delay_load = TRUE)
  pd = reticulate::import("pandas", delay_load = TRUE)
  np = reticulate::import("numpy", delay_load = TRUE)
  #mx = reticulate::import("mlx.core", delay_load = TRUE)
  if(is.null(submat)) submat = TIRTLtools::submat
  if(is.null(params)) params = TIRTLtools::params
  chunk_size = as.integer(chunk_size)
  print_chunk_size = as.integer(print_chunk_size)
  submat_py = reticulate::r_to_py(submat, convert = TRUE)
  params_py = pd$DataFrame(data = reticulate::r_to_py(params, convert = TRUE))
  tcr1_py = pd$DataFrame(data = reticulate::r_to_py(tcr1, convert = TRUE))
  tcr2_py = pd$DataFrame(data = reticulate::r_to_py(tcr2, convert = TRUE))
  res = TCRdist_gpu$TCRdist_batch(tcr1 = tcr1_py, tcr2 = tcr2_py, submat = submat_py, params_df = params_py, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_chunk_size = print_chunk_size, print_res = print_res)

  return(res)
}
