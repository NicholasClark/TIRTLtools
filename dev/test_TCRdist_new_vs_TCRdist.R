# Compare TCRdist() (R/TCRdist.R, mostly-R with python/GPU or C++ backends) against
# TCRdist_old() (original all-python implementation, kept in R/TCRdist_old.R for
# comparison) to confirm the implementations produce identical results.
#
# Run interactively, e.g.:
#   source("dev/test_TCRdist_new_vs_TCRdist.R")

pkgload::load_all(quiet = TRUE)

load_example_data(dataset = "SJTRC_minimal")
all_tcrs = get_all_tcrs(SJTRC_minimal, chain = "paired", remove_duplicates = TRUE)

## TCRdist() (new) names its returned list entries $edges_df/$nodes_df, folds tcr2
## into that single combined $nodes_df when supplied (tcr1 rows first, tcr2 rows
## numbered on after tcr1's highest tcr_index), and 1-indexes everything (R-style:
## $edges_df has "node1_idx"/"node2_idx", $nodes_df$tcr_index starts at 1).
## TCRdist_old() still uses $TCRdist_df/$tcr1/$tcr2 (a separate slot for tcr2),
## 0-indexed (Python-style): "node1_0index"/"node2_0index", tcr_index starts at 0.
## These helpers paper over all of that so the same comparison code works against
## either, by normalizing both to 1-indexed node ids.
.get_edges_df = function(res) if (!is.null(res$edges_df)) res$edges_df else res$TCRdist_df
.get_nodes_df = function(res) {
  nodes = if (!is.null(res$nodes_df)) res$nodes_df else res$tcr1
  if (!is.null(res$tcr2)) nodes = dplyr::bind_rows(nodes, res$tcr2)
  return(nodes)
}
.is_new_style = function(res) !is.null(res$edges_df)

.normalize_edges = function(res) {
  edges = as.data.frame(.get_edges_df(res))
  if (.is_new_style(res)) {
    edges = edges[, c("node1_idx", "node2_idx", "TCRdist")]
  } else {
    edges = edges[, c("node1_0index", "node2_0index", "TCRdist")]
    edges$node1_0index = edges$node1_0index + 1
    edges$node2_0index = edges$node2_0index + 1
  }
  colnames(edges) = c("node1", "node2", "TCRdist")
  edges = edges[order(edges$node1, edges$node2), ]
  rownames(edges) = NULL
  edges[] = lapply(edges, as.numeric)
  return(edges)
}

.normalize_tcr_index = function(res) {
  tcr_index = as.numeric(.get_nodes_df(res)$tcr_index)
  if (!.is_new_style(res)) tcr_index = tcr_index + 1
  return(tcr_index)
}

.compare_tcrdist_results = function(old_res, new_res, label) {
  ok = TRUE

  old_edges = .normalize_edges(old_res)
  new_edges = .normalize_edges(new_res)
  edges_match = isTRUE(all.equal(old_edges, new_edges, check.attributes = FALSE))
  if (!edges_match) ok = FALSE

  nodes_match = identical(.normalize_tcr_index(old_res), .normalize_tcr_index(new_res))
  if (!nodes_match) ok = FALSE

  status = if (ok) "PASS" else "FAIL"
  cat(sprintf(
    "[%s] %s -- edges: %d (old) vs %d (new), edges_match=%s, nodes_match=%s\n",
    status, label, nrow(old_edges), nrow(new_edges), edges_match, nodes_match
  ))
  return(ok)
}

results = c()

## Scenario 1: self-comparison, default settings, backend = "auto" (python: GPU if available, else numpy)
set.seed(1)
tcr_subset = all_tcrs[sample(nrow(all_tcrs), min(2000, nrow(all_tcrs))), ]
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, backend = "auto")
old_res = TIRTLtools:::TCRdist_old(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, print_res = FALSE, backend = "cpu")
results["self_auto_backend"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, backend='auto'")

## Scenario 2: self-comparison, backend = "cpp" (parallel C++, no python)
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400, backend = "cpp")
results["self_cpp_backend"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, backend='cpp'")

## Scenario 3: two disjoint sets of TCRs
set.seed(2)
idx1 = sample(nrow(all_tcrs), min(1200, floor(nrow(all_tcrs) / 2)))
idx2 = sample(setdiff(seq_len(nrow(all_tcrs)), idx1), min(1000, floor(nrow(all_tcrs) / 2)))
tcr1 = all_tcrs[idx1, ]
tcr2 = all_tcrs[idx2, ]
new_res = TCRdist(tcr1, tcr2, tcrdist_cutoff = 100, chunk_size = 300, backend = "auto")
old_res = TIRTLtools:::TCRdist_old(tcr1, tcr2, tcrdist_cutoff = 100, chunk_size = 300, print_res = FALSE, backend = "cpu")
results["two_sets"] = .compare_tcrdist_results(old_res, new_res, "two disjoint TCR sets")

## Scenario 4: a chunk_size that does not evenly divide n, to exercise partial chunks
new_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 333, backend = "auto")
old_res = TIRTLtools:::TCRdist_old(tcr_subset, tcrdist_cutoff = 90, chunk_size = 333, print_res = FALSE, backend = "cpu")
results["uneven_chunk_size"] = .compare_tcrdist_results(old_res, new_res, "self-comparison, uneven chunk_size")

## Scenario 5: single-chain (alpha only) TCRs
tcr_alpha_only = tcr_subset
tcr_alpha_only$cdr3b = NULL
tcr_alpha_only$vb = NULL
new_res = TCRdist(tcr_alpha_only, tcrdist_cutoff = 45, chunk_size = 400, backend = "auto")
old_res = TIRTLtools:::TCRdist_old(tcr_alpha_only, tcrdist_cutoff = 45, chunk_size = 400, print_res = FALSE, backend = "cpu")
results["single_chain_alpha"] = .compare_tcrdist_results(old_res, new_res, "single-chain (alpha only)")

cat("\n=== Summary ===\n")
print(results)
if (all(results)) {
  cat("\nAll scenarios PASSED: TCRdist() matches TCRdist_old().\n")
} else {
  stop("TCRdist() results do not match TCRdist_old() for: ", paste(names(results)[!results], collapse = ", "))
}
