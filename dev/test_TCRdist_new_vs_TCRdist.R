# Compare TCRdist_new() (R/TCRdist_new.R) against the existing TCRdist() (R/TCRdist.R)
# to confirm that the two implementations produce identical results.
#
# Run interactively, e.g.:
#   source("dev/test_TCRdist_new_vs_TCRdist.R")


pkgload::load_all(quiet = TRUE)

load_example_data(dataset = "SJTRC_minimal")
all_tcrs = get_all_tcrs(SJTRC_minimal, chain = "paired", remove_duplicates = TRUE)

.normalize_edges = function(res) {
  edges = as.data.frame(res$TCRdist_df)[, c("node1_0index", "node2_0index", "TCRdist")]
  edges = edges[order(edges$node1_0index, edges$node2_0index), ]
  rownames(edges) = NULL
  edges[] = lapply(edges, as.numeric)
  return(edges)
}

.compare_tcrdist_results = function(old_res, new_res, label) {
  ok = TRUE

  old_edges = .normalize_edges(old_res)
  new_edges = .normalize_edges(new_res)
  edges_match = isTRUE(all.equal(old_edges, new_edges, check.attributes = FALSE))
  if (!edges_match) ok = FALSE

  tcr1_match = identical(as.numeric(old_res$tcr1$tcr_index), as.numeric(new_res$tcr1$tcr_index))
  if (!tcr1_match) ok = FALSE

  tcr2_match = TRUE
  if (!is.null(old_res$tcr2) || !is.null(new_res$tcr2)) {
    tcr2_match = identical(as.numeric(old_res$tcr2$tcr_index), as.numeric(new_res$tcr2$tcr_index))
    if (!tcr2_match) ok = FALSE
  }

  status = if (ok) "PASS" else "FAIL"
  cat(sprintf(
    "[%s] %s -- edges: %d (old) vs %d (new), edges_match=%s, tcr1_match=%s, tcr2_match=%s\n",
    status, label, nrow(old_edges), nrow(new_edges), edges_match, tcr1_match, tcr2_match
  ))
  return(ok)
}

results = c()

## Scenario 1: self-comparison, default (lower triangle only)
set.seed(1)
tcr_subset = all_tcrs[sample(nrow(all_tcrs), min(2000, nrow(all_tcrs))), ]
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, print_res = FALSE, backend = "cpu")
old_res = TIRTLtools:::TCRdist_old(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, print_res = FALSE, backend = "cpu")
results["self_lower_tri"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, only_lower_tri=TRUE")

## Scenario 2: self-comparison, full (both triangles, no diagonal)
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, only_lower_tri = FALSE, print_res = FALSE, backend = "cpu")
old_res = TIRTLtools:::TCRdist_old(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, only_lower_tri = FALSE, print_res = FALSE, backend = "cpu")
results["self_full"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, only_lower_tri=FALSE")

## Scenario 3: two disjoint sets of TCRs
set.seed(2)
idx1 = sample(nrow(all_tcrs), min(1200, floor(nrow(all_tcrs) / 2)))
idx2 = sample(setdiff(seq_len(nrow(all_tcrs)), idx1), min(1000, floor(nrow(all_tcrs) / 2)))
tcr1 = all_tcrs[idx1, ]
tcr2 = all_tcrs[idx2, ]
new_res = TCRdist(tcr1, tcr2, tcrdist_cutoff = 100, chunk_size = 300, print_res = FALSE, backend = "cpu")
old_res = TIRTLtools:::TCRdist_old(tcr1, tcr2, tcrdist_cutoff = 100, chunk_size = 300, print_res = FALSE, backend = "cpu")
results["two_sets"] = .compare_tcrdist_results(old_res, new_res, "two disjoint TCR sets")

## Scenario 4: a chunk_size that does not evenly divide n, to exercise partial chunks
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 333, print_res = FALSE, backend = "cpu")
old_res = TIRTLtools:::TCRdist_old(tcr_subset, tcrdist_cutoff = 90, chunk_size = 333, print_res = FALSE, backend = "cpu")
results["uneven_chunk_size"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, uneven chunk_size")

## Scenario 5: single-chain (alpha only) TCRs
tcr_alpha_only = tcr_subset
tcr_alpha_only$cdr3b = NULL
tcr_alpha_only$vb = NULL
new_res = TCRdist(tcr_alpha_only, tcrdist_cutoff = 45, chunk_size = 400, print_res = FALSE, backend = "cpu")
old_res = TIRTLtools:::TCRdist_old(tcr_alpha_only, tcrdist_cutoff = 45, chunk_size = 400, print_res = FALSE, backend = "cpu")
results["single_chain_alpha"] = .compare_tcrdist_results(old_res, new_res, "single-chain (alpha only)")

cat("\n=== Summary ===\n")
print(results)
if (all(results)) {
  cat("\nAll scenarios PASSED: TCRdist_new() matches TCRdist().\n")
} else {
  stop("TCRdist_new() results do not match TCRdist() for: ", paste(names(results)[!results], collapse = ", "))
}
