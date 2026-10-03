# Tests for the pairing backends of run_pairing():
#   - GPU ("cupy" or "mlx"): python kernels in inst/python/pairing_new/pairing_core.py
#   - C++ ("cpp"): kernels in src/pairing_cpp.cpp
# compared against run_pairing_old() (the original python code, inst/python/pairing/) and
# against exact reference calculations in R.
#
# What "matches" means:
#   - The GPU kernels do the same float32 arithmetic as the original python code, so the GPU
#     backend must give identical results to run_pairing_old().
#   - MAD-HYPE only uses counts of wells, so the C++ backend must give identical MAD-HYPE pairs.
#   - T-SHELL correlations are computed in float32 on the GPU (abs. error up to ~1e-5) and in
#     double precision in C++. So the C++ correlations must be within 1e-5 of the GPU ones, and
#     the 3 best beta chains for an alpha chain can only differ if the beta chains are tied within
#     1e-5 (checked with exact correlations). Pairs in the final output can only differ for those
#     alpha chains, or if the adjusted p-value is close to the threshold. The p-values themselves
#     are compared where float32 is precise enough (p > 1e-50); for r close to 1, a tiny change in
#     r gives a large change in a (very small) p-value.
#
# Uses the well-level example data, which can be downloaded with:
#   TIRTLtools:::download_data("exp3_tp1_cd8.tar.gz")
# The data folder is only read from. All output is written to a temporary folder.
#
# Run interactively, e.g.:
#   source("dev/test_run_pairing_backends.R")

## mlx >= 0.31 uses reduced precision (TF32) for float32 matrix products on the GPU unless this is
## set before mlx is first used. run_pairing() sets it for the new code, so set it here for
## run_pairing_old() too, so that both compute in full float32 precision.
Sys.setenv(MLX_ENABLE_TF32 = "0")
pkgload::load_all(quiet = TRUE)
## python packages for both the old code (pandas, scipy) and the new code, before python starts
py_require(packages = get_py_deps())

data_dir = file.path(get_cache_dir(), "exp3_tp1_cd8")
if (!dir.exists(data_dir)) stop("Well-level example data not found: ", data_dir)

if (.has_nvidia_gpu()) {
  gpu_backend = "cupy"
} else if (.is_apple_silicon()) {
  gpu_backend = "mlx"
} else {
  gpu_backend = NULL ## the old code then runs on the CPU with numpy, and GPU tests are skipped
}
cat("GPU backend:", if (is.null(gpu_backend)) "none (GPU tests skipped)" else gpu_backend, "\n")

out_root = file.path(tempdir(), "test_run_pairing_backends")
unlink(out_root, recursive = TRUE)
dir.create(out_root)

r_tol = 1e-5 ## max difference between float32 (GPU) and double (C++) T-SHELL correlations
results = c()

.report = function(label, ok, detail = "") {
  status = if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", status, label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  return(ok)
}
.try_test = function(expr) {
  tryCatch(expr, error = function(e) list(ok = FALSE, detail = paste("error:", conditionMessage(e))))
}

## ===========================================================================
## Part A: kernels on simulated data
## ===========================================================================

## Simulated read fractions (chain x well) for alpha and beta chains of 'n_clones' clones, each in
## a random set of wells, plus unpaired chains. Includes a duplicated beta chain (exact ties).
.simulate_mats = function(n_wells = 48, n_clones = 150, n_unpaired = 150, seed = 1) {
  set.seed(seed)
  .chain = function(wells) { x = numeric(n_wells); x[wells] = rexp(length(wells)) * 1e-4; x }
  clone_wells = lapply(seq_len(n_clones), function(i) sample(n_wells, sample(3:25, 1)))
  alpha = t(sapply(c(clone_wells, lapply(seq_len(n_unpaired), function(i) sample(n_wells, sample(1:10, 1)))), .chain))
  beta = t(sapply(c(clone_wells, lapply(seq_len(n_unpaired), function(i) sample(n_wells, sample(1:10, 1)))), .chain))
  beta = rbind(beta, beta[5, ]) ## exact tie with beta chain 5
  rownames(alpha) = paste0("a", seq_len(nrow(alpha)))
  rownames(beta) = paste0("b", seq_len(nrow(beta)))
  colnames(alpha) = colnames(beta) = paste0("w", seq_len(n_wells))
  mats = lapply(list(alpha = alpha, beta = beta), function(m) {
    m = Matrix::Matrix(m, sparse = TRUE)
    m[Matrix::rowSums(m > 0) > 2, , drop = FALSE] ## as in run_pairing() with min_wells = 2
  })
  invisible(capture.output(mdh <- madhyper_surface(n_wells, cells = 50, alpha = 2, prior = .madhype_prior(mats))))
  return(list(mats = mats, mdh = mdh))
}

## exact MAD-HYPE pairs, computed with dense matrices
.reference_madhype = function(mats, mdh) {
  A = as.matrix(mats$alpha > 0) * 1
  B = as.matrix(mats$beta > 0) * 1
  wij = A %*% t(B)
  wa = rowSums(A)
  wb = rowSums(B)
  wi = wa - wij
  wj = t(wb - t(wij))
  thres = matrix(mdh[cbind(as.vector(wij) + 1, as.vector(wi) + 1)], nrow(wij))
  pass = which(wj < thres, arr.ind = TRUE)
  pass = pass[order(pass[, 1], pass[, 2]), , drop = FALSE]
  return(data.frame(alpha_nuc = as.integer(pass[, 1]), beta_nuc = as.integer(pass[, 2]), wij = as.numeric(wij[pass]),
                    wa = as.integer(wa[pass[, 1]]), wb = as.integer(wb[pass[, 2]])))
}

## round to float32 (exact for +, -, /: a double result rounded to float32 is the float32 result)
.f32 = function(x) readBin(writeBin(as.numeric(x), raw(), size = 4), "double", size = 4, n = length(x))

## exact T-SHELL pairs: the 3 best-correlated beta chains for each alpha chain (chains in > 2 wells),
## ties broken by the higher beta row (as for a stable ascending sort), lowest correlation first
.reference_tshell = function(mats, filter_before_top3 = FALSE) {
  wa = .n_wells_per_chain(mats$alpha)
  wb = .n_wells_per_chain(mats$beta)
  alpha_rows = which(wa > 2)
  beta_rows = which(wb > 2)
  A = as.matrix(mats$alpha[alpha_rows, , drop = FALSE])
  B = as.matrix(mats$beta[beta_rows, , drop = FALSE])
  cors = cor(t(A), t(B))
  W = (A > 0) %*% t(B > 0)
  if (filter_before_top3) {
    fwa = matrix(wa[alpha_rows], nrow(W), ncol(W))
    fwb = matrix(wb[beta_rows], nrow(W), ncol(W), byrow = TRUE)
    denom = .f32(W + (fwb - W) + (fwa - W))
    loss_sum = .f32(.f32((fwb - W) / denom) + .f32((fwa - W) / denom))
    cors[W <= 2 | loss_sum >= 0.5] = -1
  }
  n_top = min(3, ncol(cors))
  top = t(apply(cors, 1, function(r) order(r, seq_along(r))[length(r) - (n_top - 1):0]))
  if (n_top == 1) top = t(top)
  i = rep(seq_len(nrow(cors)), each = n_top)
  j = as.vector(t(top))
  alpha_nuc = alpha_rows[i]
  beta_nuc = beta_rows[j]
  return(data.frame(alpha_nuc = alpha_nuc, beta_nuc = beta_nuc, r = cors[cbind(i, j)], wij = as.numeric(W[cbind(i, j)]),
                    wa = wa[alpha_nuc], wb = wb[beta_nuc]))
}

## the original python code (as called by run_pairing_old())
.run_old_python = function(mats, mdh, filter_before_top3, chunk_size = 500L) {
  py_require(packages = get_py_deps())
  old_py = reticulate::import_from_path("pairing_all_backends", path = system.file("python/pairing/", package = "TIRTLtools"))
  invisible(capture.output(res <- old_py$pairing(
    prefix = "test", folder_out = out_root, bigmas = np_array(as.matrix(mats$alpha), dtype = "float32"),
    bigmbs = np_array(as.matrix(mats$beta), dtype = "float32"), mdh = r_to_py(mdh), write_files = FALSE,
    chunk_size = as.integer(chunk_size), filter_before_top3 = filter_before_top3)))
  if (reticulate::is_py_object(res[[1]])) res = .fix_py_to_r_df_list(res)
  return(res)
}

.same_df = function(x, y) isTRUE(all.equal(as.data.frame(x), as.data.frame(y), check.attributes = FALSE))

## Compare T-SHELL pairs from C++ (double precision) with float32 pairs (GPU). The 3 beta chains for
## an alpha chain may differ only if they are tied within r_tol (using exact correlations).
## Returns ok, detail, and the alpha chains (rows) with a different set of beta chains.
.compare_tshell_float32 = function(f32, cpp, mats) {
  problems = c()
  key = function(df) paste(df$alpha_nuc, df$beta_nuc)
  if (!identical(f32$alpha_nuc, cpp$alpha_nuc)) problems = c(problems, "different alpha chains")
  sets_f32 = split(f32$beta_nuc, f32$alpha_nuc)
  sets_cpp = split(cpp$beta_nuc, cpp$alpha_nuc)
  alphas = intersect(names(sets_f32), names(sets_cpp))
  diff_alphas = alphas[!mapply(setequal, sets_f32[alphas], sets_cpp[alphas])]
  max_gap = 0
  for (a in diff_alphas) {
    only_f32 = setdiff(sets_f32[[a]], sets_cpp[[a]])
    only_cpp = setdiff(sets_cpp[[a]], sets_f32[[a]])
    x = as.numeric(mats$alpha[as.integer(a), ])
    r_exact = function(b) sapply(b, function(j) cor(x, as.numeric(mats$beta[j, ])))
    ## C++ must have chosen the beta chains with the higher exact correlation, within r_tol
    gap = max(r_exact(only_f32)) - min(r_exact(only_cpp))
    max_gap = max(max_gap, gap)
    if (gap > r_tol) problems = c(problems, paste0("alpha ", a, ": beta chains differ by more than a near-tie (", signif(gap, 3), ")"))
  }
  shared = match(key(cpp), key(f32))
  ok_rows = !is.na(shared)
  r_diff = max(c(0, abs(cpp$r[ok_rows] - f32$r[shared[ok_rows]])))
  if (r_diff > r_tol) problems = c(problems, paste("max r difference", signif(r_diff, 3)))
  for (cl in c("wij", "wa", "wb")) {
    if (!isTRUE(all.equal(as.numeric(cpp[[cl]][ok_rows]), as.numeric(f32[[cl]][shared[ok_rows]])))) problems = c(problems, paste(cl, "differs"))
  }
  detail = sprintf("%d pairs, %d alpha chains with near-tied beta chains (max gap %.2g), max r diff %.2g",
                   nrow(cpp), length(diff_alphas), max_gap, r_diff)
  return(list(ok = length(problems) == 0, detail = paste(c(detail, problems), collapse = "; "),
              diff_alphas = as.integer(diff_alphas)))
}

sim = .simulate_mats()
cat("Simulated data:", nrow(sim$mats$alpha), "alpha chains,", nrow(sim$mats$beta), "beta chains,",
    ncol(sim$mats$alpha), "wells\n")
cpp_backend = .setup_pairing_backend("cpp", verbose = FALSE)

## Test A1: C++ MAD-HYPE is exactly the reference, for several chunk sizes
ref = .reference_madhype(sim$mats, sim$mdh)
for (chunk_size in c(1L, 7L, 1000L)) {
  cpp = .run_madhype(sim$mats, sim$mdh, cpp_backend, chunk_size, verbose = FALSE)
  results[paste0("A1_cpp_madhype_chunk", chunk_size)] = .report(
    paste("C++ MAD-HYPE = reference, chunk size", chunk_size), .same_df(ref, cpp), paste(nrow(cpp), "pairs"))
}

## Test A1b: MAD-HYPE pairs with no shared wells are called if the look-up table allows it
mdh_no_overlap = sim$mdh
mdh_no_overlap[1, ] = 4 ## call pairs with wij = 0 if wb < 4
ref = .reference_madhype(sim$mats, mdh_no_overlap)
cpp = .run_madhype(sim$mats, mdh_no_overlap, cpp_backend, 50L, verbose = FALSE)
results["A1b_cpp_madhype_no_overlap"] = .report("C++ MAD-HYPE = reference, pairs with no shared wells",
                                               .same_df(ref, cpp) && any(ref$wij == 0), paste(sum(ref$wij == 0), "pairs with wij = 0"))

## Test A2: C++ T-SHELL is the reference (exact correlations), with and without filter_before_top3
for (f in c(FALSE, TRUE)) {
  ref = .reference_tshell(sim$mats, filter_before_top3 = f)
  cpp = .run_tshell(sim$mats, cpp_backend, filter_before_top3 = f, chunk_size = 7L, verbose = FALSE)
  same_pairs = identical(paste(ref$alpha_nuc, ref$beta_nuc), paste(cpp$alpha_nuc, cpp$beta_nuc))
  same_values = isTRUE(all.equal(ref, cpp, tolerance = 1e-10, check.attributes = FALSE))
  results[paste0("A2_cpp_tshell_filter_", f)] = .report(
    paste("C++ T-SHELL = reference, filter_before_top3 =", f), same_pairs && same_values,
    paste(nrow(cpp), "pairs; same pairs in same order:", same_pairs))
}

## Test A2b: fewer than 3 beta chains
small = list(alpha = sim$mats$alpha, beta = sim$mats$beta[1:2, , drop = FALSE])
ref = .reference_tshell(small)
cpp = .run_tshell(small, cpp_backend, chunk_size = 50L, verbose = FALSE)
results["A2b_cpp_tshell_2_betas"] = .report("C++ T-SHELL = reference, 2 beta chains",
                                           isTRUE(all.equal(ref, cpp, tolerance = 1e-10, check.attributes = FALSE)))

## Test A3: GPU kernels give identical results to the original python code
if (!is.null(gpu_backend)) {
  gpu = .setup_pairing_backend(gpu_backend, verbose = FALSE)
  for (f in c(FALSE, TRUE)) {
    for (chunk_size in c(7L, 500L)) {
      cmp = .try_test({
        old = .run_old_python(sim$mats, sim$mdh, filter_before_top3 = f, chunk_size = chunk_size)
        new = .run_pairing_algorithms(sim$mats, sim$mdh, gpu, filter_before_top3 = f, chunk_size = chunk_size, verbose = FALSE)
        list(ok = .same_df(old$mdh, new$mdh) && .same_df(old$corr, new$corr),
             detail = sprintf("MAD-HYPE %s, T-SHELL %s", .same_df(old$mdh, new$mdh), .same_df(old$corr, new$corr)))
      })
      results[paste0("A3_gpu_filter_", f, "_chunk", chunk_size)] = .report(
        sprintf("GPU = old python, filter_before_top3 = %s, chunk size %d", f, chunk_size), cmp$ok, cmp$detail)
    }
  }
}

## Test A4: C++ vs the original python code (float32), on simulated data
for (f in c(FALSE, TRUE)) {
  cmp = .try_test({
    old = .run_old_python(sim$mats, sim$mdh, filter_before_top3 = f)
    cpp = .run_pairing_algorithms(sim$mats, sim$mdh, cpp_backend, filter_before_top3 = f, chunk_size = 500L, verbose = FALSE)
    ts = .compare_tshell_float32(old$corr, cpp$corr, sim$mats)
    list(ok = .same_df(old$mdh, cpp$mdh) && ts$ok, detail = paste0("MAD-HYPE identical: ", .same_df(old$mdh, cpp$mdh), "; T-SHELL: ", ts$detail))
  })
  results[paste0("A4_cpp_vs_old_filter_", f)] = .report(paste("C++ vs old python (simulated), filter_before_top3 =", f), cmp$ok, cmp$detail)
}

## Test A5: "auto" uses C++ if there is no GPU, and "cpu" is deprecated in favor of "cpp"
.with_no_gpu = function(expr) {
  ns = asNamespace("TIRTLtools")
  fns = c(".has_nvidia_gpu", ".is_apple_silicon")
  orig = mget(fns, envir = ns)
  for (fn in fns) { unlockBinding(fn, ns); assign(fn, function() FALSE, envir = ns) }
  on.exit(for (fn in fns) assign(fn, orig[[fn]], envir = ns))
  expr
}
backend_no_gpu = .with_no_gpu(.setup_pairing_backend("auto", verbose = FALSE))
results["A5_auto_no_gpu"] = .report("backend 'auto' without a GPU uses C++", backend_no_gpu$name == "cpp")
if (!is.null(gpu_backend)) {
  results["A5_auto_gpu"] = .report("backend 'auto' with a GPU uses it",
                                  .setup_pairing_backend("auto", verbose = FALSE)$name == gpu_backend)
}
old_opts = options(lifecycle_verbosity = "warning")
warned = FALSE
pb = withCallingHandlers(
  run_pairing(data_dir, file.path(out_root, "cpu_backend"), "test", wellset = get_well_subset(1:2, 1:12),
              backend = "cpu", pseudobulk_only = TRUE, verbose = FALSE),
  lifecycle_warning_deprecated = function(cnd) { warned <<- TRUE; invokeRestart("muffleWarning") })
options(old_opts)
results["A5_cpu_deprecated"] = .report("backend 'cpu' is deprecated", warned)

## ===========================================================================
## Part B: run_pairing() on the example data, compared with run_pairing_old()
## ===========================================================================

.run = function(fn, label, ...) {
  folder_out = file.path(out_root, label)
  t0 = Sys.time()
  res = suppressMessages(fn(folder_path = data_dir, folder_out = folder_out, prefix = "test", verbose = FALSE, ...))
  return(list(res = res, folder_out = folder_out, seconds = as.numeric(Sys.time() - t0, units = "secs")))
}

.read_output = function(file) {
  if (grepl("\\.rds$", file)) return(readRDS(file))
  if (grepl("\\.parquet$", file)) return(as.data.frame(nanoparquet::read_parquet(file)))
  if (grepl("\\.(tsv|tsv\\.gz|csv)$", file)) return(as.data.frame(data.table::fread(file)))
  return(readBin(file, "raw", n = file.size(file)))
}
.output_file = function(run, pattern) list.files(run$folder_out, pattern = pattern, full.names = TRUE)

## files that must be the same for all backends (everything except the T-SHELL pairs)
.compare_files = function(old, new, skip_pattern = NULL) {
  problems = c()
  files_old = sort(list.files(old$folder_out))
  files_new = sort(list.files(new$folder_out))
  if (!identical(files_old, files_new)) problems = c(problems, "different output files")
  files = intersect(files_old, files_new)
  if (!is.null(skip_pattern)) files = files[!grepl(skip_pattern, files)]
  for (f in files) {
    tol = if (grepl("_corresults\\.csv$", f)) 1e-6 else 1.5e-8 ## float32 written as text with different digits
    same = isTRUE(all.equal(.read_output(file.path(old$folder_out, f)), .read_output(file.path(new$folder_out, f)),
                            check.attributes = FALSE, tolerance = tol))
    if (!same) problems = c(problems, paste("file differs:", f))
  }
  return(problems)
}

## GPU: everything identical to run_pairing_old()
.compare_gpu_to_old = function(old, new) {
  problems = .compare_files(old, new)
  if (!.same_df(old$res, new$res)) problems = c(problems, "returned pairs differ")
  return(list(ok = length(problems) == 0, detail = paste(c(sprintf("%.0fs vs %.0fs (old)", new$seconds, old$seconds), problems), collapse = "; ")))
}

## C++: everything identical except T-SHELL pairs, which match up to float32 precision (see top of file)
.compare_cpp_to_old = function(old, new, pval_thres) {
  problems = .compare_files(old, new, skip_pattern = "_TIRTLoutput\\.tsv|_corresults\\.csv")
  ## unfiltered T-SHELL pairs (needs write_extra_files = TRUE)
  mats = list(alpha = Matrix::t(readRDS(.output_file(new, "_alpha_mat\\.rds$"))),
              beta = Matrix::t(readRDS(.output_file(new, "_beta_mat\\.rds$"))))
  ts = .compare_tshell_float32(.read_output(.output_file(old, "_corresults\\.csv$")),
                               .read_output(.output_file(new, "_corresults\\.csv$")), mats)
  if (!ts$ok) problems = c(problems, paste("unfiltered T-SHELL pairs:", ts$detail))
  near_tie_alphas = rownames(mats$alpha)[ts$diff_alphas]

  ## final pairs
  out_old = .read_output(.output_file(old, "_TIRTLoutput\\.tsv"))
  out_new = .read_output(.output_file(new, "_TIRTLoutput\\.tsv"))
  if (!.same_df(out_new, new$res)) problems = c(problems, "returned pairs differ from output file")
  mdh_old = out_old[out_old$method == "madhype", ]
  mdh_new = out_new[out_new$method == "madhype", ]
  if (!.same_df(mdh_old, mdh_new)) problems = c(problems, "MAD-HYPE pairs differ")

  ts_old = out_old[out_old$method == "tshell", ]
  ts_new = out_new[out_new$method == "tshell", ]
  ## pairs in only one output: only for near-tied alpha chains, or adjusted p-values near the threshold
  only = rbind(ts_old[!ts_old$alpha_beta %in% ts_new$alpha_beta, ], ts_new[!ts_new$alpha_beta %in% ts_old$alpha_beta, ])
  explained = only$alpha_nuc %in% near_tie_alphas | abs(log10(only$pval_adj / pval_thres)) < log10(2)
  if (!all(explained)) problems = c(problems, paste(sum(!explained), "T-SHELL pairs in only one output"))
  ## shared pairs
  shared = intersect(ts_old$alpha_beta, ts_new$alpha_beta)
  a = ts_old[match(shared, ts_old$alpha_beta), ]
  b = ts_new[match(shared, ts_new$alpha_beta), ]
  stat_cols = c("r", "ts", "pval", "pval_adj")
  if (!.same_df(a[, setdiff(colnames(a), stat_cols)], b[, setdiff(colnames(b), stat_cols)])) {
    problems = c(problems, "T-SHELL pairs differ (other than r and p-values)")
  }
  r_diff = max(c(0, abs(a$r - b$r)))
  if (r_diff > r_tol) problems = c(problems, paste("max r difference", signif(r_diff, 3)))
  precise = a$pval > 1e-50
  log_p_diff = max(c(0, abs(log10(a$pval[precise]) - log10(b$pval[precise])),
                     abs(log10(a$pval_adj[precise]) - log10(b$pval_adj[precise]))))
  if (log_p_diff > 0.05) problems = c(problems, paste("max log10 p-value difference", signif(log_p_diff, 3)))

  detail = sprintf("%.0fs vs %.0fs (old); MAD-HYPE %d pairs; T-SHELL %d/%d pairs (%d in only one, %d near-tied alpha chains), max r diff %.2g, max log10 p diff %.2g (p > 1e-50)",
                   new$seconds, old$seconds, nrow(mdh_new), nrow(ts_new), nrow(ts_old), nrow(only),
                   length(near_tie_alphas), r_diff, log_p_diff)
  return(list(ok = length(problems) == 0, detail = paste(c(detail, problems), collapse = "; ")))
}

## run the old code once, then the new code with the GPU and C++ backends, and compare
.test_backends = function(label, pval_thres, ...) {
  old = tryCatch(.run(run_pairing_old, paste0(label, "_old"), backend = if (is.null(gpu_backend)) "cpu" else gpu_backend, ...),
                 error = function(e) NULL)
  if (is.null(old)) return(.report(label, FALSE, "run_pairing_old() failed"))
  if (!is.null(gpu_backend)) {
    cmp = .try_test(.compare_gpu_to_old(old, .run(run_pairing, paste0(label, "_gpu"), backend = gpu_backend, ...)))
    results[paste0("B_", label, "_gpu")] <<- .report(paste0(label, ": ", gpu_backend, " = old"), cmp$ok, cmp$detail)
  }
  cmp = .try_test(.compare_cpp_to_old(old, .run(run_pairing, paste0(label, "_cpp"), backend = "cpp", ...), pval_thres))
  results[paste0("B_", label, "_cpp")] <<- .report(paste0(label, ": cpp vs old"), cmp$ok, cmp$detail)
}

## Test B1: default settings (192 wells pass QC, so 384-well T-SHELL settings)
.test_backends("defaults", pval_thres = 1e-10, write_extra_files = TRUE)

## Test B2: 128 wells (96-well T-SHELL settings), with all options
.test_backends("128_wells_all_options", pval_thres = 1e-3, write_extra_files = TRUE,
               wellset = get_well_subset(1:16, 1:8), filter_before_top3 = TRUE,
               select_best_madhype = TRUE, select_best_tshell = TRUE, exclude_nonfunctional = TRUE,
               gzip_output = TRUE, clone_threshold_chain = "independent")

## Test B3: user-supplied T-SHELL settings, beta clone threshold, non-default read/well filters
.test_backends("tshell_settings_filters", pval_thres = 1e-3, write_extra_files = TRUE,
               tshell_settings = get_tshell_settings("96_well"), clone_threshold_chain = "beta",
               min_reads = 1, min_wells = 3)

## Test B4: no QC (missing wells and the nearly empty well P3 are kept) and a small chunk size
.test_backends("no_qc_chunk_size", pval_thres = 1e-10, write_extra_files = TRUE, run_qc = FALSE, chunk_size = 333)

cat(sprintf("\n%d/%d tests passed\n", sum(results), length(results)))
if (!all(results)) stop("Some run_pairing() backend tests failed: ", paste(names(results)[!results], collapse = ", "))
