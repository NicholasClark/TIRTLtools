
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

.add_commas = function(x) {
 formatC(x, format="d",big.mark=",")
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

.is.DataFrame = function(data) {
  "DFrame" %in% class(data)
}

.is_df = function(data) {
  .is.DataFrame(data) || is.data.frame(data) || .is.dtplyr(data)
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

.pad_center_vec <- function(seqs, target_length) {
  sapply(seqs, function(x) .pad_center(x, target_length))
}

.prepend_newline = function(string) {
  return(paste("\n", string, sep = ""))
}

.split_string_multiline = function(string, width = 65) {
  lines = strwrap(string, width = width, simplify = TRUE)
  result = paste(lines, collapse = "\n")
  return(result)
}

.get_mode <- function(x) {
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}


.summarize_by_cdr3_nt = function(df) {
  #if(length(unique(df$targetSequences)) < nrow(df)) {
  df = df %>%
    mutate(
      allVHitsWithScore = gsub("\\*.*", "", allVHitsWithScore),
      allJHitsWithScore = gsub("\\*.*", "", allJHitsWithScore)
      ) %>%
    group_by(targetSequences) %>%
    summarize(
      aaSeqCDR3 = aaSeqCDR3[which.max(readCount)],
      allVHitsWithScore = allVHitsWithScore[which.max(readCount)],
      allJHitsWithScore = allJHitsWithScore[which.max(readCount)],
      readCount = sum(readCount),
      readFraction = sum(readFraction),
      .groups = "drop") %>%
    arrange(desc(readCount))
  #}
  return(df)
}

.get_well_from_filename = function(str_vec, well_pos=3) {
  str_vec2 = gsub(".", "_", basename(str_vec), fixed = TRUE)
  sapply(strsplit(str_vec2, split = "_"), function(x) x[[well_pos]])
}

.create_folder = function(folder_out, verbose = TRUE) {
  ## Create folder for results
  if (!dir.exists(folder_out)) {
    dir.create(folder_out, recursive = TRUE)
    msg = paste("Created folder:", folder_out)
    if(verbose) print(msg)
  } else {
    msg = paste("Folder already exists:", folder_out)
    if(verbose) print(msg)
  }
  return(invisible(NULL))
}


my_apply <- function(X, FUN, ..., parallel = FALSE, mc.cores = 2) {

  if (parallel) {
    apply_fun <- function(...) parallel::mclapply(..., mc.cores = mc.cores)
  } else {
    apply_fun <- lapply
  }

  apply_fun(X, FUN, ...)
}

calc_fisher_pval_df = function(df, n) {
  sapply(1:nrow(df), function(i) {
    m11 = df$wij[i]
    m12 = df$wi[i]
    m21 = df$wj[i]
    m22 = n - m11 - m12 - m21
    mat = matrix(c(m11, m12, m21, m22), nrow = 2)
    tmp = fisher.test(mat)
    return(tmp$p.value)
  })
}

.order_columns = function(df, cols, keep_extra = TRUE, verbose = TRUE) {
  checkmate::assert_data_frame(df)
  checkmate::assert_character(cols)
  checkmate::assert_logical(keep_extra)
  checkmate::assert_logical(verbose)

  cols_found = cols[cols %in% colnames(df)]
  cols_missing = cols[!cols %in% colnames(df)]
  extra_cols = colnames(df)[!colnames(df) %in% cols_found]
  if(verbose && length(cols_missing) >= 1L) {
    msg = paste("df does not contain columns:", paste(cols_missing, collapse = ", "))
    warning(msg)
  }
  if(keep_extra) {
    new_cols = c(cols_found, extra_cols)
  } else {
    new_cols = cols_found
  }
  if(is.data.table(df)) {
    df_out = df[, ..new_cols]
  } else {
    df_out = df[,new_cols, drop=FALSE]
  }
  return(df_out)
}

.confirm <- function(prompt, default = TRUE) {
  if (!interactive()) {
    return(default)
  }
  return(isTRUE(utils::askYesNo(prompt)))
}

pb_asset_size <- function(file, repo, tag = "latest") {
  info <- piggyback::pb_list(repo = repo, tag = tag)
  row <- info[info$file_name == file, ]
  if (nrow(row) == 0) {
    return(NA_real_)
  }
  row$size[1]  # bytes
}

format_size <- function(bytes) {
  if (is.na(bytes)) return("unknown size")
  units <- c("B", "KB", "MB", "GB")
  i <- 1
  while (bytes >= 1024 && i < length(units)) {
    bytes <- bytes / 1024
    i <- i + 1
  }
  sprintf("%.1f %s", bytes, units[i])
}

load_qs2 = function(dataset) {
  tic()
  cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
  file_short = glue("{dataset}.qs2")
  file = file.path(cache_root, "data-v1", file_short)
  if(!file.exists(file)) {
    dl = download_data(file_short)
    if(!dl) return(invisible(NULL))
  }
  msg = glue("Loading file: {file_short}...")
  message(msg)
  ts_data = qs2::qs_read(file)
  toc()
  return(ts_data)
}

get_data_dir = function(dataset, tag = "data-v1") {
  cache_path = tools::R_user_dir("TIRTLtools", which = "cache")
  path = file.path(cache_path, tag, dataset)
  if(dir.exists(path)) {
    return(path)
  } else {
    stop(glue("Directory doesn't exist: {path}"))
  }
}

# load_sjtrc_minimal = function(type = c(".qs2",".rds", ".tsv")) {
#   type = type[1]
#   if(type == ".tsv") {
#     folder = system.file("extdata/SJTRC_TIRTLseq_minimal", package = "TIRTLtools")
#     ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_")
#   } else if(type == ".rds") {
#     cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
#     file = file.path(cache_root, "data-v1", "SJTRC_minimal.rds")
#     if(!file.exists(file)) download_sjtrc_data("SJTRC_minimal.rds")
#     ts_data = readRDS(file.path(cache_root, "data-v1", "SJTRC_minimal.rds"))
#   } else {
#     tic()
#     cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
#     file = file.path(cache_root, "data-v1", "SJTRC_minimal.qs2")
#     if(!file.exists(file)) download_sjtrc_data("SJTRC_minimal.qs2")
#     file_short = basename(file)
#     msg = glue::glue("Loading file: {file_short}...")
#     message(msg)
#     ts_data = qs2::qs_read(file.path(cache_root, "data-v1", "SJTRC_minimal.qs2"))
#     toc()
#   }
#   return(ts_data)
# }

# load_sjtrc_longitudinal = function(type = c(".rds", ".tsv")) {
#   type = type[1]
#   if(type == ".tsv") {
#     folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
#     ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
#   } else if(type == ".rds") {
#     cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
#     file = file.path(cache_root, "data-v1", "SJTRC_longitudinal.rds")
#     if(!file.exists(file)) download_sjtrc_data("SJTRC_longitudinal.rds")
#     ts_data = readRDS(file)
#   } else {
#     tic()
#     cache_root <- tools::R_user_dir("TIRTLtools", which = "cache")
#     file = file.path(cache_root, "data-v1", "SJTRC_longitudinal.qs2")
#     if(!file.exists(file)) download_sjtrc_data("SJTRC_longitudinal.qs2")
#     file_short = basename(file)
#     msg = glue::glue("Loading file: {file_short}...")
#     message(msg)
#     ts_data = qs2::qs_read(file)
#     toc()
#   }
#   return(ts_data)
# }
