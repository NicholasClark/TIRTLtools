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

.tirtl_colors_gradient = function(palette = c("sea_green","sea1", "sea_multicolor"), n=20) {
  default = "sea_green"
  if(is.null(palette))  palette = default
  palette = palette[1]
  palettes = c("sea_green","sea1", "sea_multicolor")
  if(!palette %in% palettes) palette = default

  ### Color gradients
  sea1 <- colorRampPalette(c(
    "#001f23",  # Deep sea teal (very dark)
    "#014d64",  # Ocean Deep
    "#2e8c85",  # Tidepool Teal
    "#49a9a1",  # Sea Glass
    "#79d4c4",  # Seafoam
    "#b0ebe8"   # Aqua Light (brightest tint)
  ))

  sea_green <- colorRampPalette(c(
    "#003220",  # Deep algae green (very dark)
    "#1b5e4a",  # Turtle green
    "#2e8c85",  # Seaweed teal
    "#49a97a",  # Marine green
    "#79d4a2",  # Light green-seafoam
    "#b0f0c2"   # Pale seafoam green
  ))

  sea_multicolor <- colorRampPalette(c(
    "#001219",  # Abyssal Ink (deep ocean)
    "#014d64",  # Ocean Deep (dark teal)
    "#2e8c85",  # Tidepool Teal
    "#79d4c4",  # Seafoam
    "#b0ebe8",  # Aqua Light
    "#6a8e3a",  # Seagrass Green
    "#a2c523",  # Sunlit Algae
    "#f2c28b",  # Coral Sand
    "#f4b6b6",  # Shell Pink
    "#b49d79",  # Driftwood
    "#5b4e2d"   # Turtle Shell (earthy end)
  ))

  col_fun = get(palette)
  cols = do.call(col_fun, args = list(n))
  return(cols)
}

.tirtl_colors_distinct = function(palette = c("sea", "sea_alt", "tirtl"),n=Inf, verbose = FALSE) {
  default = "sea"  ## use sea as default
  if(is.null(palette)) palette = default
  palette = palette[1]
  palettes = c("sea", "sea_alt", "tirtl")
  if(!palette %in% palettes) palette = default

  ### Distinct color palettes
  sea <- c(
    "#5b4e2d",  # Turtle Shell (earthy, dark brown)
    "#49a9d9",  # Lagoon Blue (bright blue)
    "#6a8e3a",  # Seagrass Green (earthy green)
    "#fa7268",  # Reef Coral (vivid red-orange)
    "#b0ebe8",  # Aqua Light (bright cyan)
    "#395c3b",  # Kelp Forest (deep green)
    "#f2c28b",  # Coral Sand (warm beach tone)
    "#2e8c85",  # Tidepool Teal (blue-green)
    "#b24030",  # Crustacean Red (burnt red)
    "#79d4c4",  # Seafoam (pale green-cyan)
    "#78804b",  # Olive Turtle (neutral olive)
    "#014d64",  # Ocean Deep (dark teal)
    "#f4b6b6",  # Shell Pink (soft pink)
    "#a2c523",  # Sunlit Algae (bright yellow-green)
    "#cfd8d7",  # Foam Grey (cool light gray)
    "#b49d79",  # Driftwood (desaturated beige)
    "#1c5d99",  # Marine Blue (rich mid blue)
    "#e6d2ae",  # Soft Sand (light tan)
    "#001219",  # Abyssal Ink (near-black blue)
    "#f5f3e7"   # Shell White (warm off-white)
  )

  sea_alt = c(
    "#014d64",  # Ocean Deep
    "#79d4c4",  # Seafoam
    "#5b4e2d",  # Turtle Shell
    "#f2c28b",  # Coral Sand
    "#6a8e3a",  # Seagrass Green
    "#fa7268",  # Reef Coral
    "#49a9d9",  # Lagoon Blue
    "#b49d79",  # Driftwood
    "#f5f3e7",  # Shell White
    "#2e8c85",  # Tidepool Teal
    "#78804b",  # Olive Turtle
    "#1c5d99",  # Marine Blue
    "#395c3b",  # Kelp Forest
    "#e6d2ae",  # Soft Sand
    "#b0ebe8",  # Aqua Light
    "#f4b6b6",  # Shell Pink
    "#001219",  # Abyssal Ink
    "#b24030",  # Crustacean Red
    "#cfd8d7",  # Foam Grey
    "#a2c523"   # Sunlit Algae
  )

  tirtl = c(
    "#beab5c",
    "#87c529",
    "#b4ccd5",
    "#a0ad61",
    "#25d8b9",
    "#6c4f3a",
    "#346c04",
    "#e9dfae",
    "#222d4d",
    "#e2600f",
    "#4e5930",
    "#5d5b6f",
    "#64bc51",
    "#94d092",
    "#9ba91a",
    "#7e7e77",
    "#8c9cbc"
    )

  trek1 = c("#8B799CFF", "#3C999CFF", "#C86C18FF", "#E2ED50FF",
            "#13A4EBFF", "#FFF7A3FF", "#944D40FF", "#524559FF",
            "#CA480DFF", "#2F7270FF", "#F9AB3CFF",
            "#2E7BC5FF", "#BFCAFEFF", "#66FFFFFF", "#B46356FF",
            "#7A4B42FF", "#D78017FF", "#8BEAFFFF",
            "#9B5928FF", "#A1B3E2FF", "#FFE705FF")

  cols = get(palette)
  if(n>length(cols)) {
    if(verbose) warning(paste("Returning all", length(cols), "colors"))
    n=length(cols)
  }
  return(cols[1:n])
}

.get_proportion_column = function(proportion_column, is_paired) {
  if(proportion_column == "auto") {
    if(is_paired) {
      proportion_column = "wij"
    } else {
      proportion_column = "readFraction"
    }
    msg = paste("\n", "Using ", proportion_column ," for 'proportion_column'", sep = "")
    cat(msg)
  }
}

.get_value_type = function(data) {
  is_data_frame = is.data.frame(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])

  if(value_type == "auto") {
    if(is_paired) {
      value_type = "n_wells"
    } else {
      value_type = "readFraction"
    }
  }
}
