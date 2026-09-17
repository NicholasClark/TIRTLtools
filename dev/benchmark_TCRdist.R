# Benchmark TCRdist() (new, mostly-R implementation in R/TCRdist.R, with both a
# GPU/python backend and a CPU-only C++ backend) against TCRdist_old() (original
# all-python implementation, kept in R/TCRdist_old.R for comparison) at increasing
# numbers of TCRs.
#
# Run interactively, e.g.:
#   source("dev/benchmark_TCRdist.R")
#
# Adjust `sizes` and `chunk_size` below for a quicker/longer run.

pkgload::load_all(quiet = TRUE)

load_example_data(dataset = "SJTRC_longitudinal")
all_tcrs = get_all_tcrs(SJTRC_longitudinal, chain = "paired", remove_duplicates = TRUE)

sizes = c(10000, 50000, 100000)
chunk_size = 1000
tcrdist_cutoff = 90

set.seed(42)

.time_call = function(fn, ...) {
  t0 = Sys.time()
  res = fn(...)
  elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  return(list(elapsed = elapsed, res = res))
}

## Each implementation to benchmark: a label and the function call to time.
## "auto" picks a GPU backend (cupy/mlx) if available, otherwise numpy.
implementations = list(
  TCRdist_old      = function(tcr1) TIRTLtools:::TCRdist_old(tcr1 = tcr1, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_res = FALSE),
  TCRdist_auto     = function(tcr1) TCRdist(tcr1 = tcr1, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_res = FALSE, backend = "auto"),
  TCRdist_cpp      = function(tcr1) TCRdist(tcr1 = tcr1, tcrdist_cutoff = tcrdist_cutoff, chunk_size = chunk_size, print_res = FALSE, backend = "cpp")
)

benchmark_results = data.frame(
  n_tcrs_requested = integer(0),
  implementation = character(0),
  n_tcrs_after_prep = integer(0),
  n_edges = integer(0),
  elapsed_sec = numeric(0)
)

for (n in sizes) {
  n_sample = min(n, nrow(all_tcrs))
  if (n_sample < n) {
    cli::cli_alert_warning("Requested {n} TCRs but only {nrow(all_tcrs)} are available in this dataset -- sampling {n_sample} instead.")
  }
  tcr_subset = all_tcrs[sample(nrow(all_tcrs), n_sample), ]

  cli::cli_h2("n = {n} TCRs")

  for (impl_name in names(implementations)) {
    cli::cli_alert_info("Running {impl_name} ...")
    timing = .time_call(implementations[[impl_name]], tcr1 = tcr_subset)
    cli::cli_alert_success("{impl_name}: {round(timing$elapsed, 2)} sec, {nrow(timing$res$TCRdist_df)} edges")

    benchmark_results = rbind(
      benchmark_results,
      data.frame(
        n_tcrs_requested = n, implementation = impl_name,
        n_tcrs_after_prep = nrow(timing$res$tcr1),
        n_edges = nrow(timing$res$TCRdist_df),
        elapsed_sec = timing$elapsed
      )
    )
  }
}

cat("\n=== Benchmark results ===\n")
print(benchmark_results, row.names = FALSE)

wide = reshape(
  benchmark_results[, c("n_tcrs_requested", "implementation", "elapsed_sec")],
  idvar = "n_tcrs_requested", timevar = "implementation", direction = "wide"
)
colnames(wide) = gsub("^elapsed_sec\\.", "", colnames(wide))
wide$speedup_auto_vs_old = wide$TCRdist_old / wide$TCRdist_auto
wide$speedup_cpp_vs_old = wide$TCRdist_old / wide$TCRdist_cpp

cat("\n=== Speedup (TCRdist_old time / new-backend time) ===\n")
print(wide, row.names = FALSE)

out_csv = file.path("dev", "benchmark_TCRdist_results.csv")
utils::write.csv(benchmark_results, out_csv, row.names = FALSE)
cat("\nSaved raw results to", out_csv, "\n")
