# Compare run_pairing() (R/run_pairing.R, refactored) against run_pairing_old()
# (original implementation, kept in R/run_pairing_old.R for comparison) to confirm
# both produce identical results, using the python GPU backend ("cupy" or "mlx").
# Also tests the bug fixes in run_pairing() (well_pos, deprecated arguments,
# get_good_wells_new() well order check) and that pseudobulk_only skips pairing.
#
# Uses the well-level example data, which can be downloaded with:
#   TIRTLtools:::download_data("exp3_tp1_cd8.tar.gz")
# The data folder is only read from. All output is written to a temporary folder.
#
# Run interactively, e.g.:
#   source("dev/test_run_pairing.R")

pkgload::load_all(quiet = TRUE)

data_dir = file.path(get_cache_dir(), "exp3_tp1_cd8")
if (!dir.exists(data_dir)) stop("Well-level example data not found: ", data_dir)

if (.has_nvidia_gpu()) {
  gpu_backend = "cupy"
} else if (.is_apple_silicon()) {
  gpu_backend = "mlx"
} else {
  stop("No GPU found: these tests compare results using the GPU backend ('cupy' or 'mlx')")
}
cat("Using backend:", gpu_backend, "\n")

out_root = file.path(tempdir(), "test_run_pairing")
unlink(out_root, recursive = TRUE)
dir.create(out_root)

## make sure every lifecycle deprecation warning is raised (not just once per session)
old_opts = options(lifecycle_verbosity = "warning")

results = c()

.report = function(label, ok, detail = "") {
  status = if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", status, label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  return(ok)
}

## run a pairing function quietly, writing output to its own folder under out_root
.run = function(fn, label, ...) {
  folder_out = file.path(out_root, label)
  res = suppressMessages(fn(folder_out = folder_out, prefix = "test", verbose = FALSE, ...))
  return(list(res = res, folder_out = folder_out))
}

## read one output file in a form that can be compared with all.equal()
.read_output = function(file) {
  if (grepl("\\.rds$", file)) return(readRDS(file))
  if (grepl("\\.parquet$", file)) return(as.data.frame(nanoparquet::read_parquet(file)))
  if (grepl("\\.(tsv|tsv\\.gz|csv)$", file)) return(as.data.frame(data.table::fread(file)))
  return(readBin(file, "raw", n = file.size(file)))
}

## compare returned pairs and every file written to the output folders
.compare_runs = function(old, new) {
  problems = c()
  res_match = isTRUE(all.equal(as.data.frame(old$res), as.data.frame(new$res), check.attributes = FALSE))
  if (!res_match) problems = c(problems, "returned data frame differs")

  files_old = sort(list.files(old$folder_out))
  files_new = sort(list.files(new$folder_out))
  if (!identical(files_old, files_new)) {
    problems = c(problems, paste0("output files differ: old = [", paste(files_old, collapse = ", "),
                                  "], new = [", paste(files_new, collapse = ", "), "]"))
  }
  for (f in intersect(files_old, files_new)) {
    same = isTRUE(all.equal(.read_output(file.path(old$folder_out, f)),
                            .read_output(file.path(new$folder_out, f)),
                            check.attributes = FALSE))
    if (!same) problems = c(problems, paste("file differs:", f))
  }
  return(list(ok = length(problems) == 0, detail = paste(c(paste(length(files_new), "files"), problems), collapse = "; ")))
}

.test_old_vs_new = function(label, ...) {
  cmp = tryCatch({
    old = .run(run_pairing_old, paste0(label, "_old"), ...)
    new = .run(run_pairing, paste0(label, "_new"), ...)
    .compare_runs(old, new)
  }, error = function(e) list(ok = FALSE, detail = paste("error:", conditionMessage(e))))
  return(.report(label, cmp$ok, cmp$detail))
}

## run expr, capturing (and muffling) all lifecycle deprecation warnings
.capture_deprecations = function(expr) {
  msgs = c()
  val = withCallingHandlers(expr, lifecycle_warning_deprecated = function(cnd) {
    msgs <<- c(msgs, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  return(list(value = val, messages = msgs))
}
.warned_about = function(dep, arg) any(grepl(paste0("`", arg, "`"), dep$messages, fixed = TRUE))

## ===========================================================================
## Test 1: default settings. The example data is the left half of a 384-well
## plate (192 wells), so the other 192 wells in the default wellset are missing.
## ===========================================================================
results["defaults"] = .test_old_vs_new("defaults", folder_path = data_dir, backend = gpu_backend)

## ===========================================================================
## Test 2: 128 wells (< 150 wells, so "auto" picks 96-well T-SHELL settings),
## with all optional output files and the optional filtering steps.
## ===========================================================================
wells_128 = get_well_subset(1:16, 1:8)
results["128_wells_all_options"] = .test_old_vs_new(
  "128_wells_all_options", folder_path = data_dir, backend = gpu_backend,
  wellset = wells_128,
  write_extra_files = TRUE, exclude_nonfunctional = TRUE,
  select_best_madhype = TRUE, select_best_tshell = TRUE,
  gzip_output = TRUE, clone_threshold_chain = "independent"
)

## ===========================================================================
## Test 3: user-supplied T-SHELL settings, beta clone threshold, non-default
## read/well filters.
## ===========================================================================
results["tshell_settings_filters"] = .test_old_vs_new(
  "tshell_settings_filters", folder_path = data_dir, backend = gpu_backend,
  tshell_settings = get_tshell_settings("96_well"),
  clone_threshold_chain = "beta", min_reads = 1, min_wells = 3
)

## ===========================================================================
## Test 4: no QC, so the missing wells and the nearly empty well P3 are kept.
## (Note: run_qc = FALSE with min_reads = 1 fails in both versions, because well
## P3 has no beta chains with > 1 read and is dropped from the beta matrix only.)
## ===========================================================================
results["no_qc"] = .test_old_vs_new("no_qc", folder_path = data_dir, backend = gpu_backend, run_qc = FALSE)

## ===========================================================================
## Test 5: pseudobulk only. Output files should match the old function, and
## nothing should be returned.
## ===========================================================================
results["pseudobulk_only"] = .test_old_vs_new(
  "pseudobulk_only", folder_path = data_dir, backend = gpu_backend, pseudobulk_only = TRUE
)
results["pseudobulk_only_extra_files"] = .test_old_vs_new(
  "pseudobulk_only_extra_files", folder_path = data_dir, backend = gpu_backend,
  pseudobulk_only = TRUE, write_extra_files = TRUE
)

## ===========================================================================
## Test 6: pseudobulk_only should not build the read count matrices unless
## extra files are requested (big_merge_freqs2() is only used to build them).
## ===========================================================================
.count_big_merge_calls = function(...) {
  counter = new.env()
  counter$n = 0
  trace("big_merge_freqs2", where = asNamespace("TIRTLtools"), print = FALSE,
        tracer = bquote(assign("n", get("n", envir = .(counter)) + 1, envir = .(counter))))
  on.exit(untrace("big_merge_freqs2", where = asNamespace("TIRTLtools")))
  .run(run_pairing, ...)
  return(counter$n)
}
n_no_extra = .count_big_merge_calls("pseudobulk_only_count_a", folder_path = data_dir, pseudobulk_only = TRUE)
n_extra = .count_big_merge_calls("pseudobulk_only_count_b", folder_path = data_dir, pseudobulk_only = TRUE,
                                 write_extra_files = TRUE)
results["pseudobulk_only_skips_matrices"] = .report(
  "pseudobulk_only skips count matrices", n_no_extra == 0 && n_extra == 2,
  paste0("big_merge_freqs2() calls: ", n_no_extra, " without extra files, ", n_extra, " with extra files")
)

## ===========================================================================
## Test 7: well_pos is used to find the well ID in file names. Link a quarter
## plate of files into a temporary folder with an extra "X_" at the start of each
## name (so the well is in position 4), and compare with the original files.
## ===========================================================================
wells_quarter = get_well_subset(1:8, 1:12)
files_all = list.files(data_dir)
files_quarter = files_all[.get_well_from_filename(files_all, 3) %in% wells_quarter]
renamed_dir = file.path(out_root, "renamed_input")
dir.create(renamed_dir)
file.symlink(file.path(data_dir, files_quarter), file.path(renamed_dir, paste0("X_", files_quarter)))

old = .run(run_pairing_old, "well_pos_old", folder_path = data_dir, wellset = wells_quarter, backend = gpu_backend)
new = .run(run_pairing, "well_pos_new", folder_path = renamed_dir, wellset = wells_quarter, well_pos = 4,
           backend = gpu_backend)
cmp = .compare_runs(old, new)
results["well_pos"] = .report("well_pos = 4 matches renamed files", cmp$ok, cmp$detail)

## ===========================================================================
## Test 8: deprecated arguments give a warning and still work.
## ===========================================================================
## compute = FALSE is the same as pseudobulk_only = TRUE
old = .run(run_pairing_old, "compute_old", folder_path = data_dir, pseudobulk_only = TRUE)
dep = .capture_deprecations(.run(run_pairing, "compute_new", folder_path = data_dir, compute = FALSE))
cmp = .compare_runs(old, dep$value)
warned = .warned_about(dep, "compute")
results["deprecated_compute"] = .report("deprecated compute = FALSE", warned && cmp$ok,
                                        paste0("warned = ", warned, "; ", cmp$detail))

## fork and shared are ignored
for (arg in c("fork", "shared")) {
  args = list(folder_path = data_dir, pseudobulk_only = TRUE)
  args[[arg]] = TRUE
  dep = .capture_deprecations(do.call(.run, c(list(run_pairing, paste0(arg, "_new")), args)))
  results[paste0("deprecated_", arg)] = .report(paste("deprecated", arg), .warned_about(dep, arg))
}

## pval_thres_tshell and wij_thres_tshell are the same as passing tshell_settings.
## With 128 wells, "auto" would choose the 96-well settings, so these override it.
custom_tshell = list(pval_thres_tshell = 1e-6, wij_thres_tshell = 4)
old = .run(run_pairing_old, "tshell_args_old", folder_path = data_dir, backend = gpu_backend,
           wellset = wells_128, tshell_settings = custom_tshell)
dep = .capture_deprecations(
  .run(run_pairing, "tshell_args_new", folder_path = data_dir, backend = gpu_backend,
       wellset = wells_128,
       pval_thres_tshell = custom_tshell$pval_thres_tshell,
       wij_thres_tshell = custom_tshell$wij_thres_tshell)
)
cmp = .compare_runs(old, dep$value)
warned = .warned_about(dep, "pval_thres_tshell") && .warned_about(dep, "wij_thres_tshell")
results["deprecated_tshell_thresholds"] = .report(
  "deprecated pval_thres_tshell / wij_thres_tshell", warned && cmp$ok,
  paste0("warned = ", warned, "; ", cmp$detail)
)

## ===========================================================================
## Test 9: get_good_wells_new() stops with a clear error if alpha and beta wells
## are not in the same order.
## ===========================================================================
well_list = function(wells) setNames(lapply(wells, function(w) data.frame(x = 1:5)), wells)
err = tryCatch(get_good_wells_new(well_list(c("A1", "A2")), well_list(c("A2", "A1")),
                                  thres_alpha = 1, thres_beta = 1, wellset = c("A1", "A2"), verbose = FALSE),
               error = function(e) conditionMessage(e))
results["good_wells_order_check"] = .report(
  "get_good_wells_new() well order check", identical(err, "Alpha and beta wells not in the same order"),
  paste("error:", err)
)

options(old_opts)

cat(sprintf("\n%d/%d tests passed\n", sum(results), length(results)))
if (!all(results)) stop("Some run_pairing() tests failed: ", paste(names(results)[!results], collapse = ", "))
