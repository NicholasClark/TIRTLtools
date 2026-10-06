#' Find TCRalpha/beta pairs from individual well read counts
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function runs the MAD-HYPE and T-SHELL algorithms (with efficient
#' GPU implementations) on multi-well TCR-seq data, returning TCR alpha/beta chain
#' pairs predicted to come from the same clone.
#'
#' The MAD-HYPE (Multicell Analytical Deconvolution
#' for High Yield Paired-chain Evaluation) algorithm uses Bayesian
#' inference over well co-occurrence counts to call TCR alpha/beta chain
#' pairings from multicell-per-well sequencing data.
#'
#' The T-SHELL algorithm (TCRαβ Sequence Highly Efficient Linkage
#' Learning) is a heuristic for pairing TCR chains from multi-well data,
#' based on correlation of read frequency across wells.
#'
#' @details
#'
#' The function also performs QC on the wells, dropping wells with unusually small
#' numbers of clones, and calculates "pseudobulk" for TCRα and TCRβ, aggregating reads
#' over all wells. It writes three tab-separated files (.tsv), one for TCRα pseudobulk,
#' one for TCRβ pseudobulk, and one for predicted TCRαβ pairs.
#'
#' Output may include alpha or beta chains with many predicted partners. To keep only the
#' best partners by p-value for each algorithm (allowing up to two alphas for each beta chain),
#' you may run with the options `select_best_madhype = TRUE` and `select_best_tshell == TRUE`.
#' Removing nonfunctional chains with `exclude_nonfunctional = TRUE` may also help to reduce multi-pairing.
#'
#' See Pogorelyy & Kirk et al. (2026) for a full description of the T-SHELL algorithm.
#'
#' See Holec et al. (2019) for a full description of the MAD-HYPE algorithm.
#'
#'
#' @param folder_path the path of the folder with well-level data
#' @param folder_out the path of the folder to write results to. The function will create the folder if it does not exist.
#' @param prefix a prefix for the output file names
#' @param well_filter_thres wells are removed if they have fewer unique clones than: wellfilter_thres*(Avg. # of unique clones per well).
#' The default value is 0.5
#' @param min_reads minimum number of reads a chain must have in a well to be considered observed
#' (note: actual minimum is min_reads+1. Default value is 0, i.e. chain must have >= 1 read in a well)
#' Wells with no TCRalpha or no TCRbeta chains passing this filter are not used for pairing.
#' @param min_wells a chain must be observed in more than `min_wells` wells to be paired
#' (default value is 2, i.e. chain must be observed in >= 3 wells)
#' @param well_pos the position of the well ID (e.g. "B5") in the file names. For example, files
#' named "<lab>_<project>_<well_id>_TCRalpha.tsv" would use well_pos=3. (default is 3)
#' @param tshell_settings a named list with pairing settings for T-SHELL: `pval_thres_tshell` and `wij_thres_tshell`. You may pass your own list or call \code{\link{get_tshell_settings}()} with setting "auto", "384_well", or "96_well". By default, "auto" is selected and "384_well" settings will be used if there are >= 150 wells passing QC, otherwise "96_well" settings will be used.
#' - `pval_thres_tshell` is the adjusted p-value threshold for T-SHELL significance (default 1e-10 for 384-well plate, 1e-3 for 96-well plate)
#' - `wij_thres_tshell` is the threshold for the number of wells containing both chains for T-SHELL significance (default >2 wells for 384-well plate, >3 wells for 96-well plate)
#' @param wellset a vector of wells to use for the pairing
#' @param backend the computing backend for the pairing algorithms: "auto" (default) uses a GPU if one is
#' available and otherwise C++. "cupy" (NVIDIA GPU) and "mlx" (Apple silicon GPU) run on the GPU with python.
#' "cpp" runs on the CPU with parallel C++ and does not need python. ("cpu" is deprecated, use "cpp".)
#' The GPU backends compute correlations for T-SHELL in single precision (float32) and the C++ backend
#' in double precision, so T-SHELL results may differ very slightly between them.
#' @param verbose whether to print out messages (default TRUE)
#' @param write_extra_files whether to write un-necessary intermediate files (default FALSE)
#' @param filter_before_top3 whether to filter by loss fraction before extracting top 3 correlation values for T-SHELL (default FALSE)
#' @param chunk_size number of TCRalpha chains to pair at a time (default 500)
#' @param exclude_nonfunctional whether to exclude non-functional chains before pairing (default is FALSE)
#' @param select_best_madhype whether to use a secondary algorithm on the pairs from the MAD-HYPE algorithm to select
#' the best pairs for each clone (default is FALSE)
#' @param select_best_tshell whether to use a secondary algorithm on the pairs from the T-SHELL algorithm to select
#' the best pairs for each clone (default is FALSE)
#' @param gzip_output whether to compress (gzip) output files (default is FALSE)
#' @param pseudobulk_only only output pseudobulk with no pairing (default is FALSE).
#' If `write_extra_files` is also TRUE, the plate statistics and read count matrices are written too.
#' @param run_qc whether to run the quality control pipeline on the samples, dropping wells with a number of unique alpha or beta clones below a certain threshold.
#' @param clone_threshold_chain the chain to use to determine the clone threshold for QC.
#' If "alpha", clone counts for TCRalpha will be used to determine a threshold for both
#' alpha and beta chains. If "beta", TCRbeta clone counts will be used to determine
#' the threshold for both chains. If "independent", then TCRalpha clones will be used
#' to determine a QC threshold for alpha and TCRbeta clones will be used to determine
#' a separate QC threshold for beta clones. ("alpha" is default)
#' @param wellset1 `r lifecycle::badge("deprecated")` Use `wellset` instead.
#' @param compute `r lifecycle::badge("deprecated")` Use `pseudobulk_only` instead
#' (`compute = FALSE` is the same as `pseudobulk_only = TRUE`).
#' @param pval_thres_tshell `r lifecycle::badge("deprecated")` Use `tshell_settings` instead.
#' @param wij_thres_tshell `r lifecycle::badge("deprecated")` Use `tshell_settings` instead.
#' @param fork `r lifecycle::badge("deprecated")` No longer used.
#' @param shared `r lifecycle::badge("deprecated")` No longer used.
#'
#' @return
#' A data frame with the TCR-alpha/TCR-beta pairs.
#'
#' The function also writes three files to "folder_out":
#'  - A data frame ("<prefix>_pseudobulk_TRA.tsv") of pseudobulk counts for TCRalpha chains
#'  - A data frame ("<prefix>_pseudobulk_TRB.tsv") of pseudobulk counts for TCRbeta chains
#'  - A data frame ("<prefix>_TIRTLoutput.tsv") of TCR-alpha/TCR-beta pairs.
#'
#'  These files can be loaded using the \code{\link{load_tirtlseq}()} function.
#'
#'  If write_extra_files is TRUE, the function also writes sparse matrices of per-well read
#'  counts (well x chain) for TCR-alpha and beta to "<prefix>_alpha_mat.rds" and "<prefix>_beta_mat.rds".
#'  Metadata for the chains in these matrices are written to "<prefix>_alpha_meta.parquet" and
#'  "<prefix>_beta_meta.parquet" and metadata for the wells is written to "<prefix>_well_meta.parquet".
#'
#'  These files can be loaded using the \code{\link{load_well_counts_binary}()} function.
#'
#' @references
#' Pogorelyy MV, Kirk AM, Adhikari S, Minervina AA, Sundararaman B, Vegesana K,
#' Brice DC, Scott ZB, SJTRC Study Team, Thomas PG (2026). "TIRTL-seq: deep,
#' quantitative and affordable paired TCR repertoire sequencing."
#' \emph{Nature Methods}, \strong{23}, 56-64.
#' \doi{10.1038/s41592-025-02907-9}
#'
#' Holec PV, Berleant J, Bathe M, Birnbaum ME (2019). "A Bayesian framework
#' for high-throughput T cell receptor pairing." \emph{Bioinformatics},
#' \strong{35}(8), 1318-1325. \doi{10.1093/bioinformatics/bty801}
#'
#' @family pairing
#' @export
#'

run_pairing = function(
    folder_path,
    folder_out,
    prefix,
    well_filter_thres=0.5,
    min_reads=0,
    min_wells=2,
    well_pos=3,
    tshell_settings = get_tshell_settings(format = "auto"),
    wellset = get_well_subset(1:16,1:24),
    wellset1 = lifecycle::deprecated(),
    compute = lifecycle::deprecated(),
    backend = c("auto", "cupy", "mlx", "cpp"),
    pval_thres_tshell = lifecycle::deprecated(),
    wij_thres_tshell= lifecycle::deprecated(),
    verbose = TRUE,
    write_extra_files = FALSE,
    filter_before_top3 = FALSE,
    fork = lifecycle::deprecated(),
    shared = lifecycle::deprecated(),
    chunk_size = 500,
    exclude_nonfunctional = FALSE,
    select_best_madhype = FALSE,
    select_best_tshell = FALSE,
    gzip_output = FALSE,
    pseudobulk_only = FALSE,
    run_qc = TRUE,
    clone_threshold_chain = c("alpha", "beta", "independent")
){
  tictoc::tic()
  on.exit(tictoc::toc(quiet = !verbose), add = TRUE)
  clone_threshold_chain = clone_threshold_chain[1]
  checkmate::assert_choice(clone_threshold_chain, choices = c("alpha", "beta", "independent"))
  checkmate::assert_directory_exists(folder_path)

  ## 0. Arguments ===========================================================
  backend = backend[1]
  if(identical(backend, "cpu")) {
    lifecycle::deprecate_warn(
      when = "0.2.7",
      what = I('run_pairing(backend = "cpu")'),
      with = I('run_pairing(backend = "cpp")')
    )
    backend = "cpp"
  }
  checkmate::assert_choice(backend, choices = .pairing_backends)
  if (lifecycle::is_present(wellset1)) {
    lifecycle::deprecate_warn(
      when = "0.1.22",
      what = "run_pairing(wellset1=)",
      with = "run_pairing(wellset=)"
    )
    wellset = wellset1
  }
  if (lifecycle::is_present(compute)) {
    lifecycle::deprecate_warn(
      when = "0.2.2",
      what = "run_pairing(compute=)",
      with = "run_pairing(pseudobulk_only=)"
    )
    pseudobulk_only = !compute
  }
  ## user-supplied T-SHELL thresholds override 'tshell_settings' (applied after "auto" settings are chosen)
  tshell_overrides = list()
  if (lifecycle::is_present(pval_thres_tshell)) {
    lifecycle::deprecate_warn(
      when = "0.1.22",
      what = "run_pairing(pval_thres_tshell=)",
      with = "run_pairing(tshell_settings=)"
    )
    tshell_overrides$pval_thres_tshell = pval_thres_tshell
  }
  if (lifecycle::is_present(wij_thres_tshell)) {
    lifecycle::deprecate_warn(
      when = "0.1.22",
      what = "run_pairing(wij_thres_tshell=)",
      with = "run_pairing(tshell_settings=)"
    )
    tshell_overrides$wij_thres_tshell = wij_thres_tshell
  }
  if (lifecycle::is_present(fork)) {
    lifecycle::deprecate_warn(
      when = "0.1.24",
      what = "run_pairing(fork=)",
      details = "This argument is no longer used and will be ignored."
    )
  }
  if (lifecycle::is_present(shared)) {
    lifecycle::deprecate_warn(
      when = "0.1.24",
      what = "run_pairing(shared=)",
      details = "This argument is no longer used and will be ignored."
    )
  }
  if(!is.null(tshell_settings)) .validate_tshell_settings(tshell_settings)

  .create_folder(folder_out, verbose = verbose)

  ## 1. Load well data ======================================================
  well_meta = .list_well_files(folder_path, wellset, well_pos, verbose)
  counts = list(
    alpha = .read_well_files(folder_path, well_meta$file_alpha, well_meta$well, "alpha", verbose),
    beta = .read_well_files(folder_path, well_meta$file_beta, well_meta$well, "beta", verbose)
  )
  if(exclude_nonfunctional) {
    if(verbose) message("Excluding chains with stop codons or frameshifts...")
    counts = lapply(counts, function(mlist) lapply(mlist, .get_functional))
  }
  if(verbose) {
    .report_chain_totals(counts$alpha, "alpha")
    .report_chain_totals(counts$beta, "beta")
  }

  ## 2. Well QC =============================================================
  clone_thres = .get_clone_thresholds(counts, well_filter_thres, clone_threshold_chain)
  qc = NULL
  if(run_qc) {
    qc = get_good_wells_new(counts$alpha, counts$beta,
                            thres_alpha = clone_thres$alpha,
                            thres_beta = clone_thres$beta,
                            pos = well_pos,
                            wellset = wellset,
                            verbose = verbose)
  }
  if(write_extra_files) .write_plate_stats(counts, well_meta, qc, folder_out, prefix)
  if(run_qc) {
    counts$alpha = counts$alpha[qc$a]
    counts$beta = counts$beta[qc$b]
  }

  ## 3. Pseudobulk ==========================================================
  if(verbose) message("Tabulating pseudobulk counts")
  pseudobulk = lapply(counts, .make_pseudobulk)
  .write_pseudobulk(pseudobulk, folder_out, prefix, gzip_output, verbose)
  if(verbose) message("Pseudobulk done")
  if(pseudobulk_only && !write_extra_files) return(invisible(NULL))

  ## 4. Read count matrices (chain x well) ==================================
  ## both matrices have the same wells (columns), in the same order
  pairing_wells = .get_pairing_wells(counts, min_reads, verbose)
  mats = list(
    alpha = .build_count_matrix(counts$alpha, "alpha", pairing_wells, min_reads, min_wells, verbose),
    beta = .build_count_matrix(counts$beta, "beta", pairing_wells, min_reads, min_wells, verbose)
  )
  if(write_extra_files) .write_count_matrices(mats, pseudobulk, folder_out, prefix, verbose)
  if(pseudobulk_only) return(invisible(NULL))

  ## 5. Pairing =============================================================
  tshell_settings = .resolve_tshell_settings(tshell_settings, n_wells = length(counts$alpha),
                                             overrides = tshell_overrides, verbose = verbose)
  n_wells = ncol(mats$alpha)
  prior = .madhype_prior(mats)
  if(verbose) message("Pre-computing look-up table:")
  mdh = madhyper_surface(n_wells = n_wells, cells = clone_thres$alpha, alpha = 2, prior = prior)
  if(write_extra_files) write_dat(mdh, fname = file.path(folder_out, .pairing_output_file(prefix, "mdh.tsv")))

  backend = .setup_pairing_backend(backend, verbose)
  pair_res = .run_pairing_algorithms(mats, mdh, backend, filter_before_top3 = filter_before_top3,
                                     chunk_size = chunk_size, verbose = verbose)
  if(write_extra_files) .write_unfiltered_pairs(pair_res, folder_out, prefix, verbose)

  ## 6. Post-processing =====================================================
  if(verbose) message("Filtering results, adding amino acid and V segment information...")
  result = .combine_pair_results(pair_res, mats)
  if(verbose) message("Scoring unique pairs...")
  result = .score_pairs(result, n_wells = n_wells, cpw = clone_thres$alpha, prior = prior)
  result = .filter_tshell_pairs(result, tshell_settings)
  result = .annotate_pairs(result, counts)
  ### (optional) secondary algorithm to select best MAD-HYPE/T-SHELL clones
  result = .filter_pairs(result, filter_madhype = select_best_madhype, filter_tshell = select_best_tshell)

  ## 7. Write output ========================================================
  .write_pairs(result, folder_out, prefix, gzip_output, verbose)
  return(result)
}
