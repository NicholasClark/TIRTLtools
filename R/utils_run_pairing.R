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

## sparse matrix of read fractions (chain x well), keeping chains seen in more than min_wells wells
.build_count_matrix = function(mlist, chain, min_reads = 0, min_wells = 2, verbose = TRUE) {
  if(verbose) message(paste0("Merging ", chain, " clonesets..."))
  mat = big_merge_freqs2(mlist, min_reads = min_reads)
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

## run MAD-HYPE and T-SHELL in python. Returns a list with two data frames:
## "mdh" (MAD-HYPE) and "corr" (T-SHELL), with chains given by row indices of the matrices.
.run_pairing_python = function(mats, mdh, prefix, folder_out, backend, filter_before_top3,
                               chunk_size, write_files, verbose = TRUE) {
  py_require(packages = get_py_deps())
  if(verbose) message("Running pairing algorithms...")
  py_path = system.file("python/pairing/", package = "TIRTLtools")
  pairing = reticulate::import_from_path("pairing_all_backends", path = py_path, convert = TRUE, delay_load = TRUE)
  pair_res = pairing$pairing(
    prefix = prefix, folder_out = folder_out,
    bigmas = np_array(as.matrix(mats$alpha), dtype = "float32"),
    bigmbs = np_array(as.matrix(mats$beta), dtype = "float32"),
    mdh = r_to_py(mdh), backend = backend, filter_before_top3 = filter_before_top3,
    chunk_size = chunk_size, write_files = write_files
  )
  if(verbose) message("Pairing is finished.")
  ### fix for when r_to_py doesn't convert data frames
  if(reticulate::is_py_object(pair_res[[1]])) pair_res = .fix_py_to_r_df_list(pair_res)
  return(pair_res)
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
