

convert_to_summarized_experiment = function(data) {
  meta_df = data$meta
  print("binding paired data")
  paired_df = .bind_all_samples(data, "paired")
  print("binding alpha data")
  alpha_df = .bind_all_samples(data, "alpha")
  print("binding beta data")
  beta_df = .bind_all_samples(data, "beta")

  print("alpha df to SE")
  alpha_se = .df_to_SE_pseudobulk(alpha_df, meta_df)
  print("beta df to SE")
  beta_se = .df_to_SE_pseudobulk(beta_df, meta_df)
  print("paired df to SE")
  paired_se = .df_to_SE_paired(paired_df, meta_df)
  se_list <- SimpleList(
    alpha = alpha_se,
    beta = beta_se,
    paired = paired_se
  )
  return(se_list)
}

# save_TIRTL_to_HDF5 = function(se_list, dir, prefix) {
#   saveHDF5SummarizedExperiment(se_list$alpha, dir = dir, prefix = paste(prefix, "alpha", sep = "_"))
#   saveHDF5SummarizedExperiment(se_list$beta, dir = dir, prefix = paste(prefix, "beta", sep = "_"))
#   saveHDF5SummarizedExperiment(se_list$paired, dir = dir, prefix = paste(prefix, "paired", sep = "_"))
#   return(invisible(NULL))
# }
#
# load_TIRTL_HDF5 = function(dir, prefix) {
#   alpha_se = loadHDF5SummarizedExperiment(dir = dir, prefix = paste(prefix, "alpha", sep = "_"))
#   beta_se = loadHDF5SummarizedExperiment(dir = dir, prefix = paste(prefix, "beta", sep = "_"))
#   paired_se = loadHDF5SummarizedExperiment(dir = dir, prefix = paste(prefix, "paired", sep = "_"))
#   se_list <- SimpleList(
#     alpha = alpha_se,
#     beta = beta_se,
#     paired = paired_se
#   )
#   return(se_list)
# }


save_TIRTL_meta = function(se_list, prefix, chains = c("alpha", "beta", "paired"), type = c("parquet", "HDF5")) {
  type = type[1]
  sample_meta = colData(se_list[[1]])
  if(type == "parquet") {
    #saveRDS(sample_meta, file = paste(prefix,"_col_meta.rds",sep=""))
    write_parquet(sample_meta %>% as_tibble(), sink = paste(prefix,"_col_meta.parquet",sep=""))
    for(chain in chains) {
      message(paste("Saving ", chain, " metadata...", sep = ""))
      tcr_meta = rowData(se_list[[chain]]) %>% as_tibble()
      write_parquet(tcr_meta, sink = paste(prefix, "_", chain, "_row_meta.parquet", sep = ""))
    }
  } else {
    cols = colnames(sample_meta)
    file = paste(prefix,".h5",sep="")
    grp = "col_meta"
    if(!file.exists(file)) h5createFile(file)
    h5_df = h5ls(file)
    if(!paste("/", grp, sep = "") %in% h5_df$group) h5createGroup(file, grp)
    for(col in cols) {
      print(col)
      h5write(sample_meta[[col]], file, paste(grp, col, sep = "/"))
      #writeHDF5Array(sample_meta[[col]] %>% as.matrix(),  file, name = paste("metadata", "column", col, sep = "_"))
    }

    for(chain in chains) {
      print(chain)
      #cols = rowData(se_list[[chain]]) %>% colnames()
      cols = .get_metadata_cols(chain)
      grp1 = chain
      h5createGroup(file, grp1)
      grp = paste(chain, "row_meta", sep = "/")
      #h5createFile(file)
      h5_df = h5ls(file)
      print(h5_df)
      string = paste("/", grp, sep = "")
      if(!string %in% h5_df$group) h5createGroup(file, grp)
      print("test")
      for(col in cols) {
        print(col)
        vec = rowData(se_list[[chain]])[[col]]
        h5createDataset(
          file = file,
          dataset = paste(grp, col, sep = "/"),
          dims = length(vec),
          storage.mode = storage.mode(vec),
          size = max(nchar(vec)),       # max string length (important!)
          chunk = 1e5      # number of elements per chunk
        )
        h5write(vec, file, paste(grp, col, sep = "/"))
        #writeHDF5Array(rowData(se_list[[chain]])[[col]] %>% as.matrix(),  paste(prefix,".h5",sep=""), name = paste("metadata", "row", chain, col, sep = "_"))
      }
    }
  }
  return(invisible(NULL))
}

.tenx_to_dgc = function(tenx) {
  as(tenx, "dgCMatrix")
}

save_TIRTL_to_HDF5_10x = function(se_list, prefix, chains = c("alpha", "beta", "paired")) {
  for(chain in chains) {
    message(paste("Saving ", chain, " assays...", sep = ""))
    assays = assays(se_list[[chain]]) %>% names()
    for(val in assays) {
      print(val)
      writeTENxMatrix(se_list[[chain]]@assays@data[[val]], filepath = paste(prefix,".h5",sep=""),
                      group = paste(chain,val,sep="_"))
    }
  }
  save_TIRTL_meta(se_list, prefix, chains)
  return(invisible(NULL))
}

.get_metadata_cols = function(chain) {
  if(chain == "paired") {
    cols = c("alpha_nuc", "beta_nuc", "cdr3a", "va", "ja", "cdr3b", "vb", "jb")
  } else {
    cols = c("targetSequences", "v", "j", "aaSeqCDR3", "row_id")
  }
  return(cols)
}

.get_assay_names = function(chain) {
  if(chain == "paired") {
    cols = c("wi", "wj", "wij", "wa", "wb", "is_paired", "is_paired_madhype", "is_paired_tshell")
  } else {
    cols = c("readCount", "readFraction", "sem", "n_wells", "readCount_max",
             "readCount_median", "max_wells")
  }
  return(cols)
}

load_TIRTL_HDF5_10x = function(prefix, meta_type = c("parquet", "HDF5")) {
  meta_type = meta_type[1]
  file = paste(prefix, ".h5", sep = "")
  vals = .get_assay_names("single")
  a_list = lapply(vals, function(val) {
    grp = paste("alpha", val, sep = "_")
    TENxMatrix(file, group = grp)
  }) %>% set_names(vals)
  b_list = lapply(vals, function(val) {
    grp = paste("beta", val, sep = "_")
    TENxMatrix(file, group = grp)
  }) %>% set_names(vals)
  vals = .get_assay_names("paired")
  paired_list = lapply(vals, function(val) {
    grp = paste("paired", val, sep = "_")
    TENxMatrix(file, group = grp)
  }) %>% set_names(vals)
  msg = paste("Reading metadata...")
  if(meta_type == "parquet") {
    df_list = read_meta_parquet(prefix)
    ll = df_list$df
  } else {
    ll = read_meta_HDF5(prefix)
  }
  se_alpha = SummarizedExperiment::SummarizedExperiment(colData = ll$col, rowData = ll$alpha, assays = a_list)
  metadata(se_alpha)$row_data_file <- df_list$file[["alpha"]]
  se_beta = SummarizedExperiment::SummarizedExperiment(colData = ll$col, rowData = ll$beta, assays = b_list)
  metadata(se_beta)$row_data_file <- df_list$file[["beta"]]
  se_paired = SummarizedExperiment::SummarizedExperiment(colData = ll$col, rowData = ll$paired, assays = paired_list)
  metadata(se_paired)$row_data_file <- df_list$file[["paired"]]
  se_list <- S4Vectors::SimpleList(
    alpha = se_alpha,
    beta = se_beta,
    paired = se_paired
  )
  return(se_list)
}

get_row_data = function(se, backend = c("duckdb", "arrow", "duckplyr")) {
  backend = backend[1]
  file = metadata(se)$row_data_file
  if(backend == "arrow") df = arrow::open_dataset(file)
  if(backend == "duckdb") {
    con <- duckdb::dbConnect(duckdb::duckdb())
    call = paste("read_parquet(", "'", file, "'", ")", sep = "")
    #tbl <- tbl(con, "read_parquet('alpha_test.parquet')")
    df <- dplyr::tbl(con, call)
  }
  if(backend == "duckplyr") df = duckplyr::read_parquet_duckdb(file)
  return(df)
}

show_row_data = function(se, n=6, backend = c("duckdb", "arrow")) {
  df_tmp = arrow::open_dataset(metadata(se)$row_data_file)
  dims = dim(df_tmp) %>% sapply(function(x) format(x, big.mark=","))
  backend = backend[1]
  df = get_row_data(se, backend)
  txt <- format(df %>% head(n) %>% collect() )
  dim_str = paste0(dims, collapse = " × ")
  txt[1] = paste("\033[38;5;246m# A tibble: ", dim_str, "\033[39m", sep = "")
  cat(paste(txt, collapse = "\n"), "\n")
}

# read_metaRDS = function(prefix) {
#   a_file = paste(prefix, "_alpha_row_meta.rds", sep = "")
#   b_file = paste(prefix, "_beta_row_meta.rds", sep = "")
#   #p_file = paste(prefix, "_paired_row_meta.rds", sep = "")
#   col_file = paste(prefix, "_col_meta.rds", sep = "")
#   #files = c(a_file, b_file, p_file, col_file)
#   #list_names = c("alpha", "beta", "paired", "col")
#   files = c(a_file, b_file, col_file)
#   list_names = c("alpha", "beta", "col")
#   ll = lapply(files, function(x) {
#     message(paste("Loading ", x, "...", sep = ""))
#     readRDS(x)
#     }) %>% set_names(list_names)
#   return(ll)
# }

read_meta_parquet = function(prefix, cols = "row_id") {
  a_file = paste(prefix, "_alpha_row_meta.parquet", sep = "")
  b_file = paste(prefix, "_beta_row_meta.parquet", sep = "")
  p_file = paste(prefix, "_paired_row_meta.parquet", sep = "")
  col_file = paste(prefix, "_col_meta.parquet", sep = "")
  files = c(a_file, b_file, p_file, col_file)
  list_names = c("alpha", "beta", "paired", "col")
  #files = c(a_file, b_file, col_file)
  #list_names = c("alpha", "beta", "col")
  ll = lapply(1:length(files), function(i) {
    message(paste("Loading ", list_names[i], "...", sep = ""))
    if(list_names[i] %in% c("alpha", "beta")) {
      df = arrow::read_parquet(files[i], col_select = cols)
      #df = DataFrame(df)
      #df$targetSequences = DNAStringSet(df$targetSequences)
      #df$aaSeqCDR3 = BStringSet(df$aaSeqCDR3)
    } else {
      df = arrow::read_parquet(files[i])
    }
  }) %>% set_names(list_names)
  names(files) = list_names
  return(list(df = ll, files = files))
}

read_meta_HDF5 = function(prefix) {
  h5 = paste(prefix, ".h5", sep = "")
  df_alpha = read_h5_df(h5, "alpha/row_meta")
  df_beta = read_h5_df(h5, "beta/row_meta")
  df_paired = read_h5_df(h5, "paired/row_meta")
  df_col_meta = read_h5_df(h5, "col_meta")
  #df_alpha = read_pattern_h5_df(h5, "metadata_row_alpha")
  #df_beta = read_pattern_h5_df(h5, "metadata_row_beta")
  #df_paired = read_pattern_h5_df(h5, "metadata_row_paired")
  #df_col_meta = read_pattern_h5_df(h5, "metadata_column")
  ll = list(alpha = df_alpha, beta = df_beta, paired = df_paired, col = df_col_meta)
  files = rep(h5, 4) %>% set_names(c("alpha","beta","paired", "col"))
  return(ll)
}

read_h5_df = function(h5, grp) {
  cols = h5ls(h5) %>% filter(group == paste("/", grp, sep = "")) %>% extract2("name")
  ll = lapply(cols, function(col) {
    HDF5DataFrame::HDF5ColumnVector(h5, name = grp, column = col)
  })
  ddf = DelayedDataFrame::DelayedDataFrame(ll)
  colnames(ddf) = cols
  return(ddf)
}

# read_pattern_h5_df = function(h5, pattern) {
#   col_h5_names = h5ls(h5) %>% as_tibble() %>% filter(grepl(pattern, name)) %>% extract2("name")
#   col_names = gsub(paste(pattern, "_", sep = ""), "", col_h5_names)
#   df = read_h5_df(h5, col_h5_names, col_names)
#   return(df)
# }
#
# read_h5_df = function(h5, h5_names, col_names) {
#   ll = lapply(h5_names, function(x) {
#     #HDF5Array(file = h5, name = x)
#     #HDF5ColumnVector(h5, name = ?? (group), column = ?? (column name))
#   })
#   ddf = DelayedDataFrame::DelayedDataFrame(ll)
#   #ddf = do.call(get("cbind", envir = asNamespace("DelayedDataFrame")), ll)
#   colnames(ddf) = col_names
#   return(ddf)
# }

save_TIRTL_to_HDF5_sparse = function(se_list, file, chains = c("alpha", "beta", "paired"), assays = NULL) {
  h5createFile(file)
  h5createGroup(file, "alpha")
  h5createGroup(file, "beta")
  h5createGroup(file, "paired")

  for(chain in chains) {
    print(chain)
    if(is.null(assays)) assays = assays(se_list[[chain]]) %>% names()
    for(val in assays) {
      print(val)
      grp = paste(chain, val, sep = "/")
      h5createGroup(file, grp)
      A = se_list[[chain]]@assays@data[[val]]
      h5write(A@x,        file, paste(grp, "x", sep = "/"))
      h5write(A@i,        file, paste(grp, "i", sep = "/"))
      h5write(A@p,        file, paste(grp, "p", sep = "/"))
      h5write(dim(A),     file, paste(grp, "Dim", sep = "/"))
    }
  }
  return(invisible(NULL))
}

load_TIRTL_HDF5_sparse = function(file) {
  info <- h5read(file, "alpha")
  vals = names(info)
  a_list = lapply(vals, function(val) {
    new("dgCMatrix",
             x = info[[val]]$x %>% as.numeric(),
             i = info[[val]]$i %>% as.integer(),
             p = info[[val]]$p %>% as.integer(),
             Dim = info[[val]]$Dim %>% as.integer())
  }) %>% set_names(vals)

  info <- h5read(file, "beta")
  vals = names(info)
  b_list = lapply(vals, function(val) {
    new("dgCMatrix",
        x = info[[val]]$x %>% as.numeric(),
        i = info[[val]]$i %>% as.integer(),
        p = info[[val]]$p %>% as.integer(),
        Dim = info[[val]]$Dim %>% as.integer())
  }) %>% set_names(vals)
  return(list(alpha = a_list, beta = b_list))
}

# load_TIRTL_SE = function(dir, prefix) {
#   se = SummarizedExperiment::SummarizedExperiment(assays = mat_list, rowData = df_row_meta, colData = meta)
# }

.bind_all_samples = function(data, chain) {
  df_list = lapply(1:length(data$data), function(i) {
    df_tmp = data$data[[i]][[chain]]
    df_tmp$sample = names(data$data)[i]
    if(chain == "paired") df_tmp = .identify_pairing_method(df_tmp)
    return(df_tmp)
  })
  if("DFrame" %in% class(df_list[[1]])) {
    df = bindROWS(df_list[[1]], df_list[-1])
  } else {
    df = bind_rows(df_list)
  }
  return(df)
}

.identify_pairing_method = function(df) {
  df$nseq_cdr3a_v_j = paste(df$alpha_nuc, df$va, df$ja, sep = "|")
  df$nseq_cdr3b_v_j = paste(df$beta_nuc, df$vb, df$jb, sep = "|")
  df$full_tcr = paste(df$va, df$ja, df$alpha_nuc, df$beta_nuc, df$vb, df$jb, sep = "|")

  df_madhype = df %>% filter(method == "madhype")
  df_tshell = df %>% filter(method == "tshell")

  df = df[!duplicated(df$full_tcr),]

  df$is_paired = TRUE
  df$is_paired_madhype = df$full_tcr %in% df_madhype$full_tcr
  df$is_paired_tshell = df$full_tcr %in% df_tshell$full_tcr
  return(df)
}

.df_to_SE_pseudobulk = function(df, meta) {
  print("making nseq_vj")
  #df$nseq_v_j = paste(df$targetSequences, df$v, df$j, sep = "|")
  df$row_id <- vctrs::vec_group_id(df[, c("targetSequences", "v", "j")])
  #dt$col_id <- vec_group_id(dt[, c("sample")])   # or more columns if composite
  print("getting rid of duplicates")
  df_row_meta = df[!duplicated(row_id),c("targetSequences", "v", "j", "aaSeqCDR3", "row_id")]
  #df_row_meta = df[!duplicated(df$nseq_v_j), c("targetSequences", "v", "j", "aaSeqCDR3","nseq_v_j")]
  #print("matching tcrs")
  #df$row = match(df$nseq_v_j, df_row_meta$nseq_v_j)
  print("matching samples")
  df$col_id = match(df$sample, meta$sample_id)
  print("making sparse matrices")
  cols_to_2d = c("readCount", "readFraction", "sem", "n_wells", "readCount_max", "readCount_median", "max_wells")
  mat_list = lapply(cols_to_2d, function(value) {
    #Matrix::sparseMatrix(x = df[[column]], i = df$row, j = df$col ,dims = c(nrow(df_row_meta), nrow(meta)))
    sparseMatrix(
      i = df$row_id,
      j = df$col_id,
      x = df[[value]],
      dims = c(max(df$row_id), max(df$col_id))
    )
  }) %>% set_names(cols_to_2d)
  #df_row_meta = df_row_meta %>% DataFrame()
  #df_row_meta$targetSequences = Biostrings::DNAStringSet(df_row_meta$targetSequences)
  #df_row_meta$aaSeqCDR3 = Biostrings::BStringSet(df_row_meta$aaSeqCDR3)
  #df_row_meta$v = factor(df_row_meta$v)
  #df_row_meta$j = factor(df_row_meta$j)
  se = SummarizedExperiment::SummarizedExperiment(assays = mat_list, rowData = df_row_meta, colData = meta)
  return(se)
  #return(list(matrices = mat_list, row_data = df_row_meta, col_data = meta))
}


## df is a long data frame w/ pseudobulk from all samples
.df_to_triplet_df_pseudobulk = function(df, meta) {
  print("making nseq_vj")
  #df$nseq_v_j = paste(df$targetSequences, df$v, df$j, sep = "|")
  df$row_id <- vctrs::vec_group_id(df[, c("targetSequences", "v", "j")])
  #dt$col_id <- vec_group_id(dt[, c("sample")])   # or more columns if composite
  print("getting rid of duplicates")
  df_row_meta = df[!duplicated(row_id),c("targetSequences", "v", "j", "aaSeqCDR3", "row_id")]
  #df_row_meta = df[!duplicated(df$nseq_v_j), c("targetSequences", "v", "j", "aaSeqCDR3","nseq_v_j")]
  #print("matching tcrs")
  #df$row = match(df$nseq_v_j, df_row_meta$nseq_v_j)
  print("matching samples")
  df$col_id = match(df$sample, meta$sample_id)
  print("making sparse matrices")
  cols_to_2d = c("readCount", "readFraction", "sem", "n_wells", "readCount_max", "readCount_median", "max_wells")
  return(df)
}

## df is a long data frame w/ paired data from all samples (duplicates in each sample removed)
.df_to_SE_paired = function(df, meta) {
  df_row_meta = df[!duplicated(df$full_tcr),]

  df$row_id = match(df$full_tcr, df_row_meta$full_tcr)
  df$col_id = match(df$sample, meta$sample_id)
  cols_to_2d = c("wi", "wj", "wij", "wa", "wb", "is_paired", "is_paired_madhype", "is_paired_tshell")
  mat_list = lapply(cols_to_2d, function(value) {
      sparseMatrix(
        i = df$row_id,
        j = df$col_id,
        x = df[[value]],
        dims = c(max(df$row_id), max(df$col_id))
      )
  }) %>% set_names(cols_to_2d)
  se = SummarizedExperiment::SummarizedExperiment(assays = mat_list, rowData = df_row_meta, colData = meta)
  return(se)
}


.se_to_df_pseudobulk = function(se_list, chain, sample) {
  se = se_list[[chain]]
  col = "readCount"
  meta_df = se@colData
  col_idx = which(meta_df$sample_id == sample)
  mat = se@assays@data[[col]]
  keep_idx = which(mat[,col_idx] != 0)
  se_sub = se[keep_idx,]
  se_df = rowData(se_sub)
  for(x in assayNames(se_sub)) {
    se_df[[x]] = assays(se_sub)[[x]][,col_idx]
  }
  cols_order = c("targetSequences", "readCount", "v", "j", "aaSeqCDR3", "n_wells",
                 "readCount_max", "readCount_median", "sem", "readFraction", "max_wells")
  se_df = se_df[order(se_df$readCount, decreasing = TRUE),cols_order]
  return(se_df)
}
