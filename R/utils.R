### prep data for tcrdist python function
prep_for_tcrdist = function(df, params=NULL) {
  if(is.null(df)) return(df)
  if(is.null(params)) params = TIRTLtools::params
  df = .add_alleles(df) # add "*01" as allele for va and vb if necessary
  df = .filter_alleles(df, params = params) # remove alleles not found in the parameter data frame
  df = df %>% filter(nchar(cdr3a) > 5, nchar(cdr3b) > 5)
  #if(!"is_functional" %in% colnames(df))
  df = identify_non_functional_seqs(df)
  df = df %>% filter(is_functional) # remove seqs w/ stop codons (*) or frameshifts (_)
  df = as.data.frame(df) ## convert to standard data frame
  return(df)
}

tcrdist_to_matrix = function(tcr_obj) {
  n_vert = dim(tcr_obj$tcr1)[1]
  dist_df = tcr_obj$TCRdist_df %>% mutate %>%
    mutate(edge1_1index = edge1_0index+1L,
           edge2_1index = edge2_0index+1L
    )
  dist_df_to_matrix(dist_df, n_vert, 1:n_vert)
}

dist_df_to_matrix = function(dist_df, n_vertices_all, idx_keep) {
  ## temporarily switch zeros with -1
  dist_df$TCRdist_mod = ifelse(dist_df$TCRdist == 0L, -1L, dist_df$TCRdist)
  ## convert to sparse matrix
  sparse_tcrdist_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
                                            x=dist_df$TCRdist_mod, symmetric = TRUE,
                                            dims = c(n_vertices_all, n_vertices_all))
  tmp = sparse_tcrdist_mat[idx_keep, idx_keep] %>% as.matrix()
  mode(tmp) <- "integer"
  #tmp[tmp==0L] = cutoff
  tmp[tmp==0L] = NA
  tmp[tmp==-1L] = 0
  return(tmp)
}

dist_obj_to_matrix = function(obj, idx_keep) {
  tmp = dist_df_to_matrix(obj$dist_df, n_vertices_all = dim(obj$df)[1], idx_keep = idx_keep)
  return(tmp)
}

get_well_subset = function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F))
}

get_wells_from_edges = function(top_left, bottom_right, return_type = c("wells", "rows_and_columns")){
  return_type = return_type[1]
  row_start = substr(top_left, 1, 1)
  col_start = as.integer( substr(top_left, 2, nchar(top_left)) )
  row_end = substr(bottom_right, 1, 1)
  col_end = as.integer( substr(bottom_right, 2, nchar(bottom_right)) )
  check1 = col_start <= col_end
  if(!check1) stop("top_left column has to be before bottom_right column")
  cols = col_start:col_end
  row_start = match(row_start, LETTERS)
  row_end = match(row_end, LETTERS)
  check2 = row_start <= row_end
  if(!check2) stop("top_left row has to be before bottom_right row")
  rows = row_start:row_end
  row_letters = LETTERS[rows]
  if(return_type == "wells") {
    out = as.vector(unlist(sapply(row_letters,function(x)paste(x,cols,sep=""),simplify = F)))
  } else {
    out = list(rows = row_letters, columns = cols)
  }
  return(out)
}

get_all_tcrs = function(data, chain = c("paired", "alpha", "beta"), remove_duplicates = TRUE) {
  chain = chain[1]
  df_all = lapply(1:length(data$data), function(i) {
    sample_df = data$data[[i]][[chain]]
    sample_df = bind_cols(sample_df, data$meta[i,])
    return(sample_df)
  }) %>% bind_rows()
  if(chain == "paired" && remove_duplicates) df_all = remove_dupes_paired(df_all)
  return(df_all)
}

#' Remove duplicates from paired chain data.
#'
#' There will be duplicates in paired chain data when a pair is called by both
#' the T-SHELL and MAD-HYPE algorithms.
#' The input 'data' needs to be either a single data frame (paired chain) or a list
#' of data frames (paired chain)
remove_dupes_paired = function(data) {
  if(is.data.frame(data)) {
    out = data[!duplicated(data$alpha_beta),]
  } else {
    out = lapply(data, function(x) {
      x[!duplicated(x$alpha_beta),]
    })
  }
  return(out)
}

### add alleles ("*01") to va and vb if necessary (needed for TCRdist)
.add_alleles = function(df) {
  has_allele_va = grepl("\\*", df$va)
  has_allele_vb = grepl("\\*", df$vb)
  df$va_orig = df$va
  df$vb_orig = df$vb
  df = df %>%
    mutate(
      va = ifelse(has_allele_va, va, paste(va, "*01", sep = "")),
      vb = ifelse(has_allele_vb, vb, paste(vb, "*01", sep = ""))
    )
  return(df)
}

.add_newline = function(string) {
  return(paste(string, "\n", sep = ""))
}

### check that va and vb are found in the parameters for tcrdist (but don't remove them)
.check_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% mutate(va_allowed = va %in% params$feature, vb_allowed = vb %in% params$feature) %>%
    mutate(va_and_vb_allowed = va_allowed & vb_allowed)
  return(df)
}

### filter data frame so that va and vb alleles are found in the parameters for tcrdist
.filter_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% filter(va %in% params$feature, vb %in% params$feature)
  return(df)
}

.get_log_labels_neg = function(x, pseudo = 1e-6, max_val = NULL, min_val = NULL, label_zero = FALSE) {
  if(is.null(max_val)) max_val = max(x+pseudo)
  #y1 = min(x+pseudo) %>% log10() %>% floor()
  if(!is.null(min_val)) {
    y1 = min_val
  } else {
    y1 = pseudo %>% log10() %>% floor()
  }
  y2 = ifelse(max_val >= 0.1, 0, -1)
  y_breaks = y1:y2
  expr_list = sapply(y_breaks, function(x) {
    if(x %in% c(0,-1,-2)) return(as.character(10^x))
    if(x == y1 && label_zero) {
      return(bquote("0"))
    } else {
      return(bquote("10"^.(as.character(x))))
    }
  })
  labels_y <- do.call("expression", expr_list)
  return(list(brks = 10^y_breaks, labels = labels_y))
}

.get_log_labels_pos = function(x) {
  max_rank = max(x) %>% log10() %>% ceiling()
  x_breaks = 0:max_rank
  expr_list = sapply(x_breaks, function(x) {
    if(x %in% c(0,1,2)) return(as.character(10^x))
    bquote("10"^.(as.character(x)))
  })
  labels_x = do.call("expression", expr_list)
  return(list(brks = 10^x_breaks, labels = labels_x))
}

.geta<-function(dtlist){
  dtlist[grepl("TRA",names(dtlist))]
}

.getb<-function(dtlist){
  dtlist[grepl("TRB",names(dtlist))]
}

.geta_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRA",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

.getb_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRB",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

.get_labels_from_col = function(meta, label_col) {
  ### get labels for samples
  if(label_col == "Sample") {
    labels = meta[[1]]
  } else {
    labels = meta[[label_col]]
  }
  if(length(unique(labels)) != dim(meta)[1]) labels = paste(1:dim(meta)[1], labels)
  return(labels)
}

# input is a string from the "AllVHitsWithScore" column, e.g. "TRBV11-1*00(1030),TRBV11-3*00(1001)"
# output is the v-segment with the highest score, i.e. "TRBV11-1*00" in this example
.get_segments = function(segment_with_score_single) {
  v_vec = strsplit(segment_with_score_single, split = ",")[[1]]
  vs = gsub("\\(.*", "", v_vec)
  scores = gsub(".*\\(([^()]*)\\).*", "\\1", v_vec) %>% as.numeric()
  max_ind = which.max(scores)
  return(as.vector(vs[max_ind]))
}

# input is a vector, a column from a data frame like "AllVHitsWithScore"
# output is a vector of the v-segments with the highest score for each
.get_segments_all = function(v_with_score) {
  sapply(v_with_score, function(x) {
    .get_segments(x)
  }) %>% as.vector()
}

.is.dtplyr = function(data) {
  num = class(data) %>% grepl("dtplyr", .) %>% sum()
  check = num > 0
  return(check)
}

.is.list.only = function(data) {
  return( is.list(data) && !is.data.frame(data) )
}

.is.paired = function(data) {
  is_data_frame = is.data.frame(data) || .is.dtplyr(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])
  return(is_paired)
}

.na_to0<-function (x) {
  x[is.na(x)]<-0
  x
}

.pad_center <- function(seq, target_length) {
  seq = unlist(strsplit(seq, split = ""))
  seq_length <- length(seq)
  if (seq_length >= target_length) {
    return(paste(seq[1:target_length], collapse = ""))
  } else {
    total_padding <- target_length - seq_length
    first_half <- seq[1:floor(seq_length / 2)]
    second_half <- seq[(floor(seq_length / 2) + 1):seq_length]
    return(paste(c(first_half, rep("_", total_padding), second_half), collapse = "" ))
  }
}

.prepend_newline = function(string) {
  return(paste("\n", string, sep = ""))
}

.split_string_multiline = function(string, width = 65) {
  lines = strwrap(string, width = width, simplify = TRUE)
  result = paste(lines, collapse = "\n")
  return(result)
}
