pad_center <- function(seq, target_length) {
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

get_log_labels_neg = function(x) {
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

get_log_labels_pos = function(x) {
  max_rank = max(x) %>% log10() %>% ceiling()
  x_breaks = 0:max_rank
  expr_list = sapply(x_breaks, function(x) {
    if(x %in% c(0,1,2)) return(as.character(10^x))
    bquote("10"^.(as.character(x)))
  })
  labels_x = do.call("expression", expr_list)
  return(list(brks = 10^x_breaks, labels = labels_x))
}

geta<-function(dtlist){
  dtlist[grepl("TRA",names(dtlist))]
}

getb<-function(dtlist){
  dtlist[grepl("TRB",names(dtlist))]
}

geta_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRA",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

getb_sm<-function(dtlist){ #same, but smaller files!
  lapply(dtlist[grepl("TRB",names(dtlist))],function(x)x[,.(readCount,readFraction,targetSequences),])
}

get_well_subset<-function(row_range=1:16,col_range=1:24){
  unlist(sapply(LETTERS[row_range],function(x)paste(x,col_range,sep=""),simplify = F))
}

is.paired = function(data) {
  is_data_frame = is.data.frame(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])
  return(is_paired)
}

is.list.only = function(data) {
  return( is.list(data) && !is.data.frame(data) )
}

get_labels_from_col = function(meta, label_col) {
  ### get labels for samples
  if(label_col == "Sample") {
    labels = meta[[1]]
  } else {
    labels = meta[[label_col]]
  }
  if(length(unique(labels)) != dim(meta)[1]) labels = paste(1:dim(meta)[1], labels)
  return(labels)
}

split_string_multiline = function(string, width = 65) {
  lines = strwrap(string, width = width, simplify = TRUE)
  result = paste(lines, collapse = "\n")
  return(result)
}

add_newline = function(string) {
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


get_type_column = function(type_column, is_paired) {
  if(type_column == "auto") {
    if(is_paired) type_column = "alpha_beta"
    if(!is_paired) type_column = "targetSequences"
    msg = paste("\n", "Using ", type_column ," for 'type_column'", sep = "")
    cat(msg)
  }
  return(type_column)
}

get_proportion_column = function(proportion_column, is_paired, is_annotated) {
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
add_alleles = function(df) {
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
filter_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% filter(va %in% params$feature, vb %in% params$feature)
  return(df)
}

### check that va and vb are found in the parameters for tcrdist (but don't remove them)
check_alleles = function(df, params=NULL) {
  if(is.null(params)) params = TIRTLtools::params
  df = df %>% mutate(va_allowed = va %in% params$feature, vb_allowed = vb %in% params$feature) %>%
    mutate(va_and_vb_allowed = va_allowed & vb_allowed)
  return(df)
}

### prep data for tcrdist python function
prep_for_tcrdist = function(df, params=NULL) {
  if(is.null(df)) return(df)
  if(is.null(params)) params = TIRTLtools::params
  df = add_alleles(df) # add "*01" as allele for va and vb if necessary
  df = filter_alleles(df, params = params) # remove alleles not found in the parameter data frame
  if(!"is_functional" %in% colnames(df)) df = identify_non_functional_seqs(df)
  df = df %>% filter(is_functional) # remove seqs w/ stop codons (*) or frameshifts (_)
  df = as.data.frame(df) ## convert to standard data frame
  return(df)
}
