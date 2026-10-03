## helper functions for run_pairing() ------------------------------------------
## Most helpers work on one chain at a time, and run_pairing() applies them to a
## list with one item per chain: list(alpha = ..., beta = ...)

## "auto" T-SHELL settings use the 384-well settings if at least this many wells pass QC
.tshell_auto_min_wells = 150

.chain_names = c(alpha = "TCRalpha", beta = "TCRbeta")

## output file name: "<prefix>_<suffix>", with ".gz" added if gzip is TRUE
.pairing_output_file = function(prefix, suffix, gzip = FALSE) {
  file = paste0(prefix, "_", suffix)
  if(gzip) file = paste0(file, ".gz")
  return(file)
}

## arguments ------------------------------------------------------------------

.validate_tshell_settings = function(tshell_settings) {
  checkmate::assert_subset(c("pval_thres_tshell", "wij_thres_tshell"), names(tshell_settings))
  checkmate::assert_number(tshell_settings$pval_thres_tshell, lower = 0, upper = 1)
  checkmate::assert_int(tshell_settings$wij_thres_tshell, lower = 0, upper = 383)
  return(invisible(tshell_settings))
}

## choose T-SHELL settings for "auto" (NULL) based on the number of wells passing QC,
## then apply any thresholds passed with the deprecated run_pairing() arguments
.resolve_tshell_settings = function(tshell_settings, n_wells, overrides = list(), verbose = TRUE) {
  if(is.null(tshell_settings)) {
    format = if(n_wells >= .tshell_auto_min_wells) "384_well" else "96_well"
    if(verbose) {
      message("Using 'auto' T-SHELL settings, selecting 'pval_thres_tshell' and 'wij_thres_tshell' based on the number of wells passing QC")
      message(paste0(n_wells, " wells passing: Using '", format, "' T-SHELL settings"))
    }
    tshell_settings = get_tshell_settings(format)
  }
  tshell_settings = utils::modifyList(tshell_settings, overrides)
  .validate_tshell_settings(tshell_settings)
  return(tshell_settings)
}

## loading well data -----------------------------------------------------------

## data frame with one row per well in 'wellset' (sorted by row and column) and
## columns 'file_alpha' and 'file_beta' with the file names (NA if missing)
.list_well_files = function(folder_path, wellset, well_pos = 3, verbose = TRUE) {
  .files_by_well = function(pattern) {
    files = list.files(path = folder_path, full.names = FALSE, pattern = pattern)
    tibble(well = .get_well_from_filename(files, well_pos), file = files) %>%
      filter(well %in% wellset)
  }
  files = list(alpha = .files_by_well(".*_TRA\\.tsv"), beta = .files_by_well(".*_TRB\\.tsv"))
  for(chain in names(files)) {
    missing_wells = wellset[!wellset %in% files[[chain]]$well]
    if(verbose && length(missing_wells) > 0) {
      message(paste0(length(missing_wells), " wells missing .tsv file (", .chain_names[[chain]], "): ",
                     paste0(missing_wells, collapse = ", ")))
    }
  }

  well_meta = tibble(well = unique(wellset)) %>%
    left_join(rename(files$alpha, file_alpha = file), by = "well") %>%
    left_join(rename(files$beta, file_beta = file), by = "well") %>%
    arrange(substr(well, 1, 1), as.integer(substr(well, 2, 3)))
  return(well_meta)
}

## read clone files into a list named by well. Wells without a file get an empty table.
.read_well_files = function(folder_path, files, wells, chain, verbose = TRUE) {
  n_files = sum(!is.na(files))
  n_loaded = cumsum(!is.na(files)) ## number of files loaded after reading file i
  if(verbose) message(paste0("Loading clone files (", .chain_names[[chain]], ") for ", n_files, " wells..."))
  mlist = lapply(seq_along(files), function(i) {
    if(is.na(files[i])) return(tibble())
    df = fread(file.path(folder_path, files[i]))
    if(verbose && (n_loaded[i] %% 25 == 0 || n_loaded[i] == n_files)) {
      message(paste(n_loaded[i], "of", n_files, "files loaded"))
    }
    return(df)
  })
  names(mlist) = wells
  return(mlist)
}

.report_chain_totals = function(mlist, chain) {
  all_wells = bind_rows(mlist)
  message(paste0("Total number of ", .chain_names[[chain]], " reads: ", scales::scientific(sum(all_wells$readCount))))
  message(paste0("Total number of unique ", .chain_names[[chain]], " chains: ", length(unique(all_wells$targetSequences))))
  return(invisible(NULL))
}

## well QC ---------------------------------------------------------------------

## QC thresholds: well_filter_thres * (avg. number of unique clones in non-empty wells)
.get_clone_thresholds = function(counts, well_filter_thres, clone_threshold_chain) {
  thres = lapply(counts, function(mlist) {
    n_clones = sapply(mlist, nrow)
    round(well_filter_thres * mean(n_clones[n_clones != 0]))
  })
  if(clone_threshold_chain == "alpha") thres$beta = thres$alpha ## use alpha threshold for both
  if(clone_threshold_chain == "beta") thres$alpha = thres$beta ## use beta threshold for both
  return(thres)
}

## read counts and number of clones in each well (before QC), with QC result if QC was run
.write_plate_stats = function(counts, well_meta, qc, folder_out, prefix) {
  .sum_reads = function(mlist) sapply(mlist, function(x) if(nrow(x) > 0) sum(x$readCount) else 0)
  stats_a = tibble(a_names = well_meta$file_alpha, well = well_meta$well,
                   a_sum_counts = .sum_reads(counts$alpha), a_rows = sapply(counts$alpha, nrow))
  stats_b = tibble(b_names = well_meta$file_beta, well = well_meta$well,
                   b_sum_counts = .sum_reads(counts$beta), b_rows = sapply(counts$beta, nrow))
  if(!is.null(qc)) {
    stats_a$qc_pass_a = qc$a
    stats_b$qc_pass_b = qc$b
  }
  plate_stats = full_join(stats_a, stats_b, by = "well")
  fwrite(plate_stats, file.path(folder_out, .pairing_output_file(prefix, "plate_stats.tsv")), sep = "\t")
  return(invisible(NULL))
}

## pseudobulk ------------------------------------------------------------------

## sum reads for each chain over all wells
.make_pseudobulk = function(mlist) {
  pseudobulk = .combineTCR(rbindlist(mlist, idcol = "file"))
  pseudobulk$max_wells = length(mlist)
  return(pseudobulk[order(-readCount), ])
}

.write_pseudobulk = function(pseudobulk, folder_out, prefix, gzip_output = FALSE, verbose = TRUE) {
  suffixes = c(alpha = "pseudobulk_TRA.tsv", beta = "pseudobulk_TRB.tsv")
  for(chain in names(pseudobulk)) {
    file = .pairing_output_file(prefix, suffixes[[chain]], gzip_output)
    if(verbose) message(paste0("Writing ", .chain_names[[chain]], " pseudobulk file... ", file))
    fwrite(pseudobulk[[chain]], file.path(folder_out, file), sep = "\t")
  }
  return(invisible(NULL))
}

## read count matrices ---------------------------------------------------------

## wells used for pairing: wells with at least one chain passing min_reads for both TCRalpha
## and TCRbeta, in the order of the well lists. A well with no chains for only one of the
## two is removed for both, so that the columns of the alpha and beta matrices match.
.get_pairing_wells = function(counts, min_reads = 0, verbose = TRUE) {
  if(!identical(names(counts$alpha), names(counts$beta))) stop("Alpha and beta wells not in the same order")
  has_chains = lapply(counts, function(mlist) {
    vapply(mlist, function(x) nrow(x) > 0 && any(x$readCount > min_reads), logical(1))
  })
  keep = has_chains$alpha & has_chains$beta
  if(verbose) {
    ## wells with no data for both chains (e.g. missing files) are not reported again
    has_data = sapply(counts$alpha, nrow) > 0 | sapply(counts$beta, nrow) > 0
    for(chain in names(has_chains)) {
      removed = names(keep)[has_data & !has_chains[[chain]]]
      if(length(removed) > 0) {
        message(paste0(length(removed), " wells with no ", .chain_names[[chain]], " chains with more than ",
                       min_reads, " reads removed from pairing (for both chains): ", paste0(removed, collapse = ", ")))
      }
    }
  }
  return(names(keep)[keep])
}

## sparse matrix of read fractions (chain x well) with columns 'wells', keeping chains seen in more than min_wells wells
.build_count_matrix = function(mlist, chain, wells, min_reads = 0, min_wells = 2, verbose = TRUE) {
  if(verbose) message(paste0("Merging ", chain, " clonesets..."))
  mat = big_merge_freqs2(mlist, min_reads = min_reads, wells = wells)
  if(verbose) message(paste("Done! Unique", chain, "clones after filtering:", nrow(mat)))
  mat = mat[rowSums(mat > 0) > min_wells, , drop = FALSE]
  if(verbose) message(paste0("Unique ", chain, " clones in more than ", min_wells, " wells: ", nrow(mat)))
  return(mat)
}

## write read count matrices (well x chain) with metadata for chains and wells
.write_count_matrices = function(mats, pseudobulk, folder_out, prefix, verbose = TRUE) {
  if(verbose) message("Writing TCRalpha/beta read count matrices (well x chain) and single-chain metadata...")
  if(!identical(colnames(mats$alpha), colnames(mats$beta))) {
    stop("Error creating alpha/beta read count matrices: wells are not in the same order")
  }
  for(chain in names(mats)) {
    file_mat = .pairing_output_file(prefix, paste0(chain, "_mat.rds"))
    file_meta = .pairing_output_file(prefix, paste0(chain, "_meta.parquet"))
    chain_meta = pseudobulk[[chain]][match(rownames(mats[[chain]]), pseudobulk[[chain]]$targetSequences), ]
    if(verbose) message(paste0("Writing ", .chain_names[[chain]], " read count matrix... ", file_mat))
    saveRDS(t(mats[[chain]]), file = file.path(folder_out, file_mat)) ## columns are clones
    if(verbose) message(paste0("Writing chain metadata for ", .chain_names[[chain]], " matrix... ", file_meta))
    nanoparquet::write_parquet(chain_meta, file.path(folder_out, file_meta))
  }
  file_well_meta = .pairing_output_file(prefix, "well_meta.parquet")
  well_meta = data.frame(idx = seq_len(ncol(mats$alpha)), well = colnames(mats$alpha))
  if(verbose) message(paste("Writing well metadata for read count matrices...", file_well_meta))
  nanoparquet::write_parquet(well_meta, file.path(folder_out, file_well_meta))
  return(invisible(NULL))
}

## pairing ---------------------------------------------------------------------

## prior probability that a given alpha and beta chain are a true pair
.madhype_prior = function(mats) {
  1/(as.numeric(nrow(mats$alpha)) * as.numeric(nrow(mats$beta)))^0.5
}

## pairing backends ------------------------------------------------------------
## The MAD-HYPE and T-SHELL computations run either on a GPU, with the python kernels in
## inst/python/pairing_new/pairing_core.py ("cupy" for NVIDIA, "mlx" for Apple silicon),
## or on the CPU with the C++ kernels in src/pairing_cpp.cpp ("cpp"). Everything else
## (inputs, chunking, outputs) is done in R and is the same for all backends.

.pairing_backends = c("auto", "cupy", "mlx", "cpp")

## Returns a list with 'name' ("cupy", "mlx", or "cpp") and 'kernels' (the python module
## for a GPU backend). "auto" uses a GPU if one is found and its python library can be
## loaded, and otherwise C++.
.setup_pairing_backend = function(backend = "auto", verbose = TRUE) {
  name = backend
  if(backend == "auto") {
    name = if(.has_nvidia_gpu()) "cupy" else if(.is_apple_silicon()) "mlx" else "cpp"
  }
  kernels = NULL
  if(name == "mlx" && !nzchar(Sys.getenv("MLX_ENABLE_TF32"))) {
    ## mlx >= 0.31 uses TF32 for float32 matrix products on the GPU by default, which makes the
    ## T-SHELL correlations accurate to only ~1e-3. Turn it off (mlx reads this once, when first used).
    Sys.setenv(MLX_ENABLE_TF32 = "0")
  }
  if(name != "cpp") {
    kernels = tryCatch({
      py_require(packages = .get_py_deps_new())
      kernels = reticulate::import_from_path("pairing_core", path = system.file("python/pairing_new/", package = "TIRTLtools"),
                                             convert = TRUE)
      kernels$check_backend(name)
      kernels
    }, error = function(e) {
      if(backend != "auto") stop("Could not load the '", name, "' backend: ", conditionMessage(e))
      if(verbose) message("Could not load the '", name, "' backend, using 'cpp' instead: ", conditionMessage(e))
      NULL
    })
    if(is.null(kernels)) name = "cpp"
  }
  if(verbose) message("Using '", name, "' backend for pairing")
  return(list(name = name, kernels = kernels))
}

## number of wells with each chain (rows of a chain x well matrix)
.n_wells_per_chain = function(mat) as.integer(Matrix::rowSums(mat > 0))

## number of wells shared by alpha chains alpha_rows and beta chains beta_rows (one value per pair)
.count_shared_wells = function(alpha, beta, alpha_rows, beta_rows) {
  shared = (alpha[alpha_rows, , drop = FALSE] > 0) & (beta[beta_rows, , drop = FALSE] > 0)
  return(as.integer(Matrix::rowSums(shared)))
}

## mean and centered L2 norm of each chain's read fractions across wells (for correlations)
.chain_stats = function(mat) {
  mat = Matrix::drop0(mat)
  n_wells = ncol(mat)
  chain_mean = Matrix::rowSums(mat) / n_wells
  centered_sq = mat
  centered_sq@x = (mat@x - chain_mean[mat@i + 1])^2 ## nonzero wells; zero wells are added below
  n_zero = n_wells - .n_wells_per_chain(mat)
  chain_norm = sqrt(Matrix::rowSums(centered_sq) + n_zero * chain_mean^2)
  return(list(mean = chain_mean, norm = chain_norm))
}

## run a pairing kernel on chunks of alpha chains (rows 1..n_alpha). kernel(rows) returns a
## list with "i" (0-based position in 'rows') and "j" (0-based beta chain) for each pair, plus
## other values. Returns a data.table with all pairs, with i and j as 1-based row numbers.
.run_kernel_in_chunks = function(kernel, n_alpha, chunk_size, label, verbose = TRUE) {
  if(n_alpha == 0) return(data.table())
  starts = seq.int(1L, n_alpha, by = chunk_size)
  milestones = .get_progress_milestones(length(starts), 10)
  if(verbose) message(paste0("Running ", label, " (", length(starts), " chunks)..."))
  res = lapply(seq_along(starts), function(k) {
    rows = starts[k]:min(starts[k] + chunk_size - 1L, n_alpha)
    pairs = lapply(kernel(rows), as.vector)
    pairs$i = rows[pairs$i + 1L]
    pairs$j = pairs$j + 1L
    if(verbose && k %in% milestones) message(paste0(label, ": ", round(100 * k / length(starts)), "% done"))
    as.data.table(pairs)
  })
  return(rbindlist(res))
}

## MAD-HYPE: all alpha/beta pairs that pass the MAD-HYPE look-up table 'mdh'. Returns a data frame
## with alpha_nuc, beta_nuc (row numbers in mats$alpha, mats$beta), wij (wells with both chains),
## and wa, wb (wells with the alpha, beta chain), sorted by alpha_nuc and beta_nuc.
.run_madhype = function(mats, mdh, backend, chunk_size = 500, verbose = TRUE) {
  wa = .n_wells_per_chain(mats$alpha)
  wb = .n_wells_per_chain(mats$beta)
  storage.mode(mdh) = "integer"
  if(backend$name == "cpp") {
    alpha_t = Matrix::drop0(Matrix::t(mats$alpha)) ## wells x alpha chains
    beta = Matrix::drop0(mats$beta)
    beta_by_wb = order(wb) - 1L ## for pairs with no shared wells
    kernel = function(rows) {
      madhype_chunk_cpp(alpha_t, beta, wa, wb, mdh, beta_by_wb, start = rows[1] - 1L, end = rows[length(rows)])
    }
  } else {
    gpu = backend$kernels$MadHype(beta_bin = np_array(as.matrix(mats$beta > 0), dtype = "float32"),
                                  b_total = np_array(wb, dtype = "int32"),
                                  mdh = np_array(mdh, dtype = "int32"), backend = backend$name)
    kernel = function(rows) {
      gpu$chunk(alpha_bin = np_array(as.matrix(mats$alpha[rows, , drop = FALSE] > 0), dtype = "float32"),
                a_total = np_array(wa[rows], dtype = "int32"))
    }
  }
  pairs = .run_kernel_in_chunks(kernel, nrow(mats$alpha), chunk_size, "MAD-HYPE", verbose)
  if(verbose) message(paste("MAD-HYPE pairs:", nrow(pairs)))
  return(data.frame(alpha_nuc = as.integer(pairs$i), beta_nuc = as.integer(pairs$j), wij = as.numeric(pairs$wij),
                    wa = wa[pairs$i], wb = wb[pairs$j]))
}

## T-SHELL: the 3 beta chains with the highest correlation of read fractions across wells for each
## alpha chain, using chains seen in more than 2 wells. Returns a data frame with alpha_nuc, beta_nuc
## (row numbers in mats$alpha, mats$beta), r (correlation), wij, wa, and wb, sorted by alpha_nuc
## and, for each alpha chain, from lowest to highest correlation.
.run_tshell = function(mats, backend, filter_before_top3 = FALSE, chunk_size = 500, verbose = TRUE) {
  wa = .n_wells_per_chain(mats$alpha)
  wb = .n_wells_per_chain(mats$beta)
  alpha_rows = which(wa > 2)
  beta_rows = which(wb > 2)
  alpha = mats$alpha[alpha_rows, , drop = FALSE]
  beta = mats$beta[beta_rows, , drop = FALSE]
  if(backend$name == "cpp") {
    alpha_t = Matrix::drop0(Matrix::t(alpha)) ## wells x alpha chains
    beta = Matrix::drop0(beta)
    alpha_stats = .chain_stats(alpha)
    beta_stats = .chain_stats(beta)
    ## for pairs with no shared wells, the correlation is highest for the beta chains with the
    ## lowest mean/norm. Ties: the beta chain with the higher row number first (same as the GPU).
    beta_order = order(beta_stats$mean / beta_stats$norm, -seq_along(beta_rows)) - 1L
    kernel = function(rows) {
      tshell_chunk_cpp(alpha_t, beta, alpha_stats$mean, alpha_stats$norm, beta_stats$mean, beta_stats$norm,
                       wa[alpha_rows], wb[beta_rows], beta_order, filter_before_top3,
                       start = rows[1] - 1L, end = rows[length(rows)])
    }
  } else {
    gpu = backend$kernels$TShell(beta = np_array(as.matrix(beta), dtype = "float32"),
                                 b_total = np_array(wb[beta_rows], dtype = "int32"),
                                 filter_before_top3 = filter_before_top3, backend = backend$name)
    kernel = function(rows) {
      gpu$chunk(alpha = np_array(as.matrix(alpha[rows, , drop = FALSE]), dtype = "float32"),
                a_total = np_array(wa[alpha_rows][rows], dtype = "int32"))
    }
  }
  pairs = .run_kernel_in_chunks(kernel, length(alpha_rows), chunk_size, "T-SHELL", verbose)
  alpha_nuc = alpha_rows[pairs$i]
  beta_nuc = beta_rows[pairs$j]
  return(data.frame(alpha_nuc = alpha_nuc, beta_nuc = beta_nuc, r = as.numeric(pairs$r),
                    wij = as.numeric(.count_shared_wells(mats$alpha, mats$beta, alpha_nuc, beta_nuc)),
                    wa = wa[alpha_nuc], wb = wb[beta_nuc]))
}

## run MAD-HYPE and T-SHELL. Returns a list with two data frames: "mdh" (MAD-HYPE) and "corr" (T-SHELL)
.run_pairing_algorithms = function(mats, mdh, backend, filter_before_top3 = FALSE, chunk_size = 500, verbose = TRUE) {
  ## the alpha and beta matrices must have the same wells (columns)
  if(!identical(colnames(mats$alpha), colnames(mats$beta))) {
    stop("Error creating alpha/beta read count matrices: wells are not in the same order")
  }
  chunk_size = as.integer(chunk_size)
  pair_res = list(
    mdh = .run_madhype(mats, mdh, backend, chunk_size, verbose),
    corr = .run_tshell(mats, backend, filter_before_top3, chunk_size, verbose)
  )
  if(verbose) message("Pairing is finished.")
  return(pair_res)
}

## write the MAD-HYPE and T-SHELL pairs before filtering and post-processing
.write_unfiltered_pairs = function(pair_res, folder_out, prefix, verbose = TRUE) {
  files = c(mdh = "madhyperesults.csv", corr = "corresults.csv")
  for(method in names(files)) {
    file = .pairing_output_file(prefix, files[[method]])
    if(verbose) message(paste("Writing unfiltered pairs...", file))
    fwrite(pair_res[[method]], file.path(folder_out, file))
  }
  return(invisible(NULL))
}

## post-processing -------------------------------------------------------------

## combine MAD-HYPE and T-SHELL pairs into one table, with chain sequences, T-SHELL
## p-values, and the fraction of wells where only one of the two chains was seen
.combine_pair_results = function(pair_res, mats) {
  n_wells = ncol(mats$alpha)
  .add_sequences = function(df) {
    df %>%
      mutate(alpha_nuc_seq = rownames(mats$alpha)[alpha_nuc],
             beta_nuc_seq = rownames(mats$beta)[beta_nuc]) %>%
      mutate(alpha_nuc = alpha_nuc_seq, beta_nuc = beta_nuc_seq,
             alpha_beta = paste0(alpha_nuc_seq, "_", beta_nuc_seq))
  }
  madhype = pair_res$mdh %>%
    .add_sequences() %>%
    mutate(method = "madhype") %>%
    as.data.table()
  ## T-SHELL results -- top 3 betas for each alpha
  tshell = pair_res$corr %>%
    .add_sequences() %>%
    mutate(ts = r*sqrt((n_wells - 2) / (1 - r^2)),
           pval = 2*pt(-abs(ts), n_wells - 2)) %>%
    as.data.table()
  tshell[, pval_adj := pval/sort(pval)[3], alpha_nuc]
  tshell[, method := "tshell"]

  result = rbind(madhype, tshell, fill = TRUE)
  ## wa/wb: wells with alpha/beta chain, wij: wells with both
  result[, loss_a_frac := (wb-wij)/(wij+(wb-wij)+(wa-wij))]
  result[, loss_b_frac := (wa-wij)/(wij+(wb-wij)+(wa-wij))]
  result[, wi := wa-wij]
  result[, wj := wb-wij]
  return(result)
}

## add MAD-HYPE score (log10 odds of a true pair) to all pairs. Each unique combination
## of well counts is scored once and merged back into the pairs.
.score_pairs = function(result, n_wells, cpw, prior) {
  well_counts = unique(result[, .(wi, wj, wij)])
  well_counts$score = vapply(seq_len(nrow(well_counts)), function(i) {
    log10(estimate_pair_prob(wi = well_counts$wi[i], wj = well_counts$wj[i], w_ij = well_counts$wij[i],
                             w_tot = n_wells, cpw = cpw, alpha = 2, prior = prior))
  }, numeric(1))
  return(merge(result, well_counts, by = c("wi", "wj", "wij"), all.x = TRUE))
}

## keep all MAD-HYPE pairs and T-SHELL pairs passing the significance thresholds
.filter_tshell_pairs = function(result, tshell_settings) {
  result[method == "madhype" |
           (method == "tshell" &
              wij > tshell_settings$wij_thres_tshell &
              pval_adj < tshell_settings$pval_thres_tshell &
              (loss_a_frac + loss_b_frac) < 0.5), ]
}

## add amino acid CDR3 and most common V and J segments for alpha and beta chains
.annotate_pairs = function(result, counts) {
  annot_a = add_VJ_aa(result$alpha_nuc, rbindlist(counts$alpha))
  result$cdr3a = annot_a$cdr3aa
  result$va = annot_a$v
  result$ja = annot_a$j

  annot_b = add_VJ_aa(result$beta_nuc, rbindlist(counts$beta))
  result$cdr3b = annot_b$cdr3aa
  result$vb = annot_b$v
  result$jb = annot_b$j
  return(result)
}

.write_pairs = function(result, folder_out, prefix, gzip_output = FALSE, verbose = TRUE) {
  file = .pairing_output_file(prefix, "TIRTLoutput.tsv", gzip_output)
  if(verbose) message(paste("Writing TCRalpha/beta pairs...", file))
  fwrite(result, file.path(folder_out, file), sep = "\t")
  if(verbose) {
    message("All pairing is finished!")
    message(paste("Number of clones paired by MAD-HYPE algorithm:", sum(result$method == "madhype")))
    message(paste("Number of clones paired by T-SHELL algorithm:", sum(result$method == "tshell")))
    message(paste("Total number of unique TCRalpha/beta pairs:", length(unique(result$alpha_beta))))
  }
  return(invisible(NULL))
}
