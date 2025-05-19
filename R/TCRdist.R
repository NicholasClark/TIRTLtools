## TCRdist main function

# tcr_all = arrow::read_parquet("../covid_data_analysis/data/vdj_and_covid_032425.parquet")
# tcr1 = tcr_all[1:20000,]
# tcr2 = tcr_all[1:20000,]
# tmp1 = TCRdist(tcr1 = tcr1, only_lower_tri = FALSE)
# tmp2 = TCRdist(tcr1 = tcr1, only_lower_tri = TRUE)
# tmp3 = TCRdist(tcr1 = tcr1, tcrdist_cutoff=120, only_lower_tri = TRUE)
# tmp2 = TCRdist(tcr1 = tcr1, tcr2 = tcr2, tcrdist_cutoff=120, only_lower_tri = F)

TCRdist = function(tcr1=NULL, tcr2=NULL, params = NULL, submat = NULL, tcrdist_cutoff=90, chunk_size=1000, print_chunk_size=1000, print_res = TRUE, only_lower_tri = TRUE) {
  ### prep R input
  # tcr1 = add_alleles(tcr1)
  # tcr1 = filter_alleles(tcr1, params = params)
  # tcr1 = as.data.frame(tcr1)
  # if(!is.null(tcr2)) tcr2 = as.data.frame(tcr2)
  tcr1 = prep_for_tcrdist(tcr1, params = params)
  tcr2 = prep_for_tcrdist(tcr2, params = params)
  chunk_size = as.integer(chunk_size)
  print_chunk_size = as.integer(print_chunk_size)
  ### load python packages/scripts
  TCRdist_gpu = reticulate::import_from_path("TCRdist_gpu", path = system.file("python", package = "TIRTLtools"), convert = TRUE, delay_load = TRUE)
  pd = reticulate::import("pandas", delay_load = TRUE)
  np = reticulate::import("numpy", delay_load = TRUE)
  ### load substitution matrix and dataframe of parameters
  if(is.null(submat)) submat = TIRTLtools::submat
  if(is.null(params)) params = TIRTLtools::params
  ### convert R objects to python objects
  submat_py = reticulate::r_to_py(submat, convert = TRUE)
  params_py = pd$DataFrame(data = reticulate::r_to_py(params, convert = TRUE))
  tcr1_py = pd$DataFrame(data = reticulate::r_to_py(tcr1, convert = TRUE))
  if(!is.null(tcr2)) {
    tcr2_py = pd$DataFrame(data = reticulate::r_to_py(tcr2, convert = TRUE))
  } else {
    tcr2_py = reticulate::r_to_py(tcr2, convert = TRUE)
  }
  ### call python TCRdist_batch function
  res = TCRdist_gpu$TCRdist_batch(tcr1 = tcr1_py, tcr2 = tcr2_py, submat = submat_py, params_df = params_py, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_chunk_size = print_chunk_size, print_res = print_res, only_lower_tri = only_lower_tri)
  return(res)
}
