#' Parallelized C++ implementation of TCRdist (no GPU required)
#'
#' @family repertoire_analysis
#'
TCRdist_cpp = function(tcr1, tcr2=NULL, tcrdist_cutoff = 90) {

  submat = TIRTLtools::submat
  params = TIRTLtools::params
  params_list = (params$value + 1L) %>% as.list() %>% set_names(params$feature)

  tcr2_equals_tcr1 = FALSE
  if(is.null(tcr2) | identical(tcr1,tcr2)) tcr2_equals_tcr1 = TRUE
  if(is.null(tcr2)) tcr2 = tcr1

  tcr1 = prep_for_tcrdist(tcr1)
  tcr2 = prep_for_tcrdist(tcr2)

  n_tcrs1 = dim(tcr1)[1]
  n_tcrs2 = dim(tcr2)[1]

  trim_len = 29

  tcr1_cdr3a = .pad_center_vec(tcr1$cdr3a, trim_len) %>% strsplit("")
  tcr1_cdr3b = .pad_center_vec(tcr1$cdr3b, trim_len) %>% strsplit("")
  tcr2_cdr3a = .pad_center_vec(tcr2$cdr3a, trim_len) %>% strsplit("")
  tcr2_cdr3b = .pad_center_vec(tcr2$cdr3b, trim_len) %>% strsplit("")

  tcr1_cdr3a_num =  sapply(1:length(tcr1_cdr3a), function(i) params_list[tcr1_cdr3a[[i]]] %>% unlist()) %>% t()
  tcr1_cdr3b_num =  sapply(1:length(tcr1_cdr3b), function(i) params_list[tcr1_cdr3b[[i]]] %>% unlist()) %>% t()

  tcr2_cdr3a_num =  sapply(1:length(tcr2_cdr3a), function(i) params_list[tcr2_cdr3a[[i]]] %>% unlist()) %>% t()
  tcr2_cdr3b_num =  sapply(1:length(tcr2_cdr3b), function(i) params_list[tcr2_cdr3b[[i]]] %>% unlist()) %>% t()

  tcr1_segments = cbind(tcr1$vb, tcr1$va)
  tcr2_segments = cbind(tcr2$vb, tcr2$va)

  tcr1_segments_num = apply(tcr1_segments, 1, function(x) params_list[x] %>% unlist()) %>% t()
  tcr2_segments_num = apply(tcr2_segments, 1, function(x) params_list[x] %>% unlist()) %>% t()

  all_mat1 = cbind(tcr1_cdr3a_num, tcr1_cdr3b_num, tcr1_segments_num)
  all_mat2 = cbind(tcr2_cdr3a_num, tcr2_cdr3b_num, tcr2_segments_num)

  tcrdist_mat <- tcrdist_parallel(all_mat1, all_mat2, submat, tcr2_equals_tcr1, tcrdist_cutoff)
  mode(tcrdist_mat) <- "integer"


  return(tcrdist_mat)
}


# TCRdist_cpp = function(tcr1, tcr2=NULL) {
#
#   submat = TIRTLtools::submat
#   params = TIRTLtools::params
#   params_list = (params$value + 1L) %>% as.list() %>% set_names(params$feature)
#
#   tcr2_equals_tcr1 = FALSE
#   if(is.null(tcr2) | identical(tcr1,tcr2)) tcr2_equals_tcr1 = TRUE
#   if(is.null(tcr2)) tcr2 = tcr1
#
#   tcr1 = prep_for_tcrdist(tcr1)
#   tcr2 = prep_for_tcrdist(tcr2)
#
#   n_tcrs1 = dim(tcr1)[1]
#   n_tcrs2 = dim(tcr2)[1]
#
#   trim_len = 29
#
#   tcr1_cdr3a = .pad_center_vec(tcr1$cdr3a, trim_len) %>% strsplit("")
#   tcr1_cdr3b = .pad_center_vec(tcr1$cdr3b, trim_len) %>% strsplit("")
#   tcr2_cdr3a = .pad_center_vec(tcr2$cdr3a, trim_len) %>% strsplit("")
#   tcr2_cdr3b = .pad_center_vec(tcr2$cdr3b, trim_len) %>% strsplit("")
#
#   tcr1_cdr3a_num =  sapply(1:length(tcr1_cdr3a), function(i) params_list[tcr1_cdr3a[[i]]] %>% unlist()) %>% t()
#   tcr1_cdr3b_num =  sapply(1:length(tcr1_cdr3b), function(i) params_list[tcr1_cdr3b[[i]]] %>% unlist()) %>% t()
#
#   tcr2_cdr3a_num =  sapply(1:length(tcr2_cdr3a), function(i) params_list[tcr2_cdr3a[[i]]] %>% unlist()) %>% t()
#   tcr2_cdr3b_num =  sapply(1:length(tcr2_cdr3b), function(i) params_list[tcr2_cdr3b[[i]]] %>% unlist()) %>% t()
#
#   tcr1_segments = cbind(tcr1$vb, tcr1$va)
#   tcr2_segments = cbind(tcr2$vb, tcr2$va)
#
#   tcr1_segments_num = apply(tcr1_segments, 1, function(x) params_list[x] %>% unlist()) %>% t()
#   tcr2_segments_num = apply(tcr2_segments, 1, function(x) params_list[x] %>% unlist()) %>% t()
#
#   all_mat1 = cbind(tcr1_cdr3a_num, tcr1_cdr3b_num, tcr1_segments_num)
#   all_mat2 = cbind(tcr2_cdr3a_num, tcr2_cdr3b_num, tcr2_segments_num)
#
#   tcrdist_mat <- tcrdist_cpp(all_mat1, all_mat2, submat, tcr2_equals_tcr1)
#   mode(tcrdist_mat) <- "integer"
#
#   return(tcrdist_mat)
# }
