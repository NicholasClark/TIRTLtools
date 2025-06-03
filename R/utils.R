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

.get_log_labels_neg = function(x) {
  y1 = min(x) %>% log10() %>% floor()
  y2 = -1
  y_breaks = y1:y2
  expr_list = sapply(y_breaks, function(x) {
    if(x %in% c(0,-1,-2)) return(as.character(10^x))
    bquote("10"^.(as.character(x)))
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

.get_well_subset<-function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F))
}

.is.dtplyr = function(data) {
  num = class(data) %>% grepl("dtplyr", .) %>% sum()
  check = num > 0
  return(check)
}

.is.paired = function(data) {
  is_data_frame = is.data.frame(data) || .is.dtplyr(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])
  return(is_paired)
}

.is.list.only = function(data) {
  return( is.list(data) && !is.data.frame(data) )
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

.split_string_multiline = function(string, width = 65) {
  lines = strwrap(string, width = width, simplify = TRUE)
  result = paste(lines, collapse = "\n")
  return(result)
}

.add_newline = function(string) {
  return(paste(string, "\n", sep = ""))
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


.get_type_column = function(type_column, is_paired) {
  if(type_column == "auto") {
    if(is_paired) type_column = "alpha_beta"
    if(!is_paired) type_column = "targetSequences"
    msg = paste("\n", "Using ", type_column ," for 'type_column'", sep = "")
    cat(msg)
  }
  return(type_column)
}

.get_proportion_column = function(proportion_column, is_paired, is_annotated) {
  if(proportion_column == "auto") {
    if(is_paired) {
      if(is.null(is_annotated) || is_annotated) {
        proportion_column = "beta_readFraction"
      } else {
        proportion_column = "wij"
      }
    } else {
      proportion_column = "readFraction"
    }
    msg = paste("\n", "Using ", proportion_column ," for 'proportion_column'", sep = "")
    cat(msg)
  }
  return(proportion_column)
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

### filter data frame so that va and vb alleles are found in the parameters for tcrdist
.filter_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% filter(va %in% params$feature, vb %in% params$feature)
  return(df)
}

### check that va and vb are found in the parameters for tcrdist (but don't remove them)
.check_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% mutate(va_allowed = va %in% params$feature, vb_allowed = vb %in% params$feature) %>%
    mutate(va_and_vb_allowed = va_allowed & vb_allowed)
  return(df)
}

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

.get_colors_25 = function() {
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  return(c25)
}

.get_colors_12 = function() {
  c12 =  c(
    "dodgerblue2", "#E31A1C",
    "green4","deeppink1",
    "#6A3D9A", "#FF7F00",
    "black", "gold1",
    "palegreen2","darkorange4",
    "orchid1", "darkturquoise")
  return(c12)
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


get_dataset_info = function(dataset = "minimal") {
  ds_options = c("minimal")
  dataset = dataset[1]
  if(!dataset %in% ds_options) stop("dataset not available")
  if(dataset == "minimal") {
    files = gh("GET /repos/:owner/:repo/contents/:path",
               owner = "NicholasClark",
               repo = "TIRTLtools_data",
               path = "SJTRC_minimal_dataset")

    file_info = lapply(files, function(x) list(name = x$name, download_url = x$download_url))
    file_info = do.call(rbind, lapply(file_info, as.data.frame))
    return(file_info)
  }
}

download_minimal_dataset = function(directory = "SJTRC_TIRTLseq_minimal") {
  file_info = get_dataset_info()
  dir.create(directory, showWarnings = FALSE)

  for (i in 1:dim(file_info)[1]) {
    url = file_info$download_url[i]
    file_name = file_info$name[i]
    msg = paste("Downloading file:", file_name, " to directory: ", directory) %>% .add_newline()
    cat(msg)
    dest = file.path(directory, basename(url))
    download.file(url, destfile = dest, mode = "wb")
  }
  return(invisible(NULL))
}
