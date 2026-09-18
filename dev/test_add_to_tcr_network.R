# Tests for add_to_tcr_network() (R/add_to_tcr_network.R).
#
# Run interactively, e.g.:
#   source("dev/test_add_to_tcr_network.R")

pkgload::load_all(quiet = TRUE)

load_example_data(dataset = "SJTRC_minimal")
all_tcrs = get_all_tcrs(SJTRC_minimal, chain = "paired", remove_duplicates = TRUE)

results = c()

.report = function(label, ok, detail = "") {
  status = if (ok) "PASS" else "FAIL"
  cat(sprintf("[%s] %s%s\n", status, label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
  return(ok)
}

.normalize_edges = function(edges_df) {
  e = as.data.frame(edges_df)[, c("node1_idx", "node2_idx", "TCRdist")]
  e = e[order(e$node1_idx, e$node2_idx), ]
  rownames(e) = NULL
  e[] = lapply(e, as.numeric)
  return(e)
}

## ===========================================================================
## Test 1: splitting a set of TCRs, building a network on the first part, then
## add_to_tcr_network()-ing the rest, should give the same result as running
## TCRdist() once on the full set.
## ===========================================================================
set.seed(1)
tcr_subset = all_tcrs[sample(nrow(all_tcrs), min(1500, nrow(all_tcrs))), ]
n1 = floor(nrow(tcr_subset) / 2)
part1 = tcr_subset[seq_len(n1), ]
part2 = tcr_subset[(n1 + 1):nrow(tcr_subset), ]

full_res = TCRdist(tcr_subset, tcrdist_cutoff = 90, chunk_size = 400)
part1_res = TCRdist(part1, tcrdist_cutoff = 90, chunk_size = 400)
combined_res = add_to_tcr_network(part1_res$edges_df, part1_res$nodes_df, part2,
                                   tcrdist_cutoff = 90, chunk_size = 400)

ok = TRUE
ok = .report("same number of nodes", nrow(combined_res$nodes_df) == nrow(full_res$nodes_df),
             sprintf("%d vs %d", nrow(combined_res$nodes_df), nrow(full_res$nodes_df))) && ok
ok = .report("same number of edges", nrow(combined_res$edges_df) == nrow(full_res$edges_df),
             sprintf("%d vs %d", nrow(combined_res$edges_df), nrow(full_res$edges_df))) && ok
ok = .report("tcr_index is a clean 1..n sequence in both",
             identical(combined_res$nodes_df$tcr_index, seq_len(nrow(combined_res$nodes_df))) &&
               identical(full_res$nodes_df$tcr_index, seq_len(nrow(full_res$nodes_df)))) && ok
ok = .report("node metadata (cdr3a/cdr3b/va/vb) identical, in the same order",
             isTRUE(all.equal(
               as.data.frame(combined_res$nodes_df[, c("tcr_index", "cdr3a", "va", "cdr3b", "vb")]),
               as.data.frame(full_res$nodes_df[, c("tcr_index", "cdr3a", "va", "cdr3b", "vb")])
             ))) && ok
ok = .report("edges (node1_idx, node2_idx, TCRdist) identical",
             isTRUE(all.equal(.normalize_edges(combined_res$edges_df), .normalize_edges(full_res$edges_df),
                               check.attributes = FALSE))) && ok
ok = .report("no duplicate edges in combined result",
             !any(duplicated(paste(combined_res$edges_df$node1_idx, combined_res$edges_df$node2_idx)))) && ok

results["split_and_rejoin_matches_full_run"] = ok

## ===========================================================================
## Test 2: the old edges_df/nodes_df rows are left completely untouched by
## add_to_tcr_network() (only appended to).
## ===========================================================================
old_nodes_untouched = identical(
  as.data.frame(combined_res$nodes_df[seq_len(nrow(part1_res$nodes_df)), ]),
  as.data.frame(part1_res$nodes_df)
)
old_edges_untouched = identical(
  as.data.frame(combined_res$edges_df[seq_len(nrow(part1_res$edges_df)), ]),
  as.data.frame(part1_res$edges_df)
)
results["old_rows_untouched"] = .report("previous edges_df/nodes_df rows unchanged, only appended to",
                                         old_nodes_untouched && old_edges_untouched)

## ===========================================================================
## Test 3: chaining multiple add_to_tcr_network() calls keeps the network
## internally consistent (clean indexing, no duplicate edges), and it stays
## usable by TCRdist_to_igraph()/TCRdist_to_sparse_matrix().
## ===========================================================================
set.seed(2)
idx = sample(nrow(all_tcrs), min(3000, nrow(all_tcrs)))
batch1 = all_tcrs[idx[1:1000], ]
batch2 = all_tcrs[idx[1001:2000], ]
batch3 = all_tcrs[idx[2001:3000], ]

chained_res = TCRdist(batch1, tcrdist_cutoff = 90)
chained_res = add_to_tcr_network(chained_res$edges_df, chained_res$nodes_df, batch2,
                                  tcrdist_cutoff = 90)
chained_res = add_to_tcr_network(chained_res$edges_df, chained_res$nodes_df, batch3,
                                  tcrdist_cutoff = 90)

once_res = TCRdist(all_tcrs[idx[1:3000], ], tcrdist_cutoff = 90)

all.equal(once_res$edges_df, chained_res$edges_df)
all.equal(once_res$nodes_df, chained_res$nodes_df)

ok = TRUE
ok = .report("chained: tcr_index is a clean 1..n sequence",
             identical(chained_res$nodes_df$tcr_index, seq_len(nrow(chained_res$nodes_df)))) && ok
ok = .report("chained: no duplicate edges",
             !any(duplicated(paste(chained_res$edges_df$node1_idx, chained_res$edges_df$node2_idx)))) && ok

gr = TCRdist_to_igraph(chained_res$edges_df, chained_res$nodes_df)
adj_mat = TCRdist_to_sparse_matrix(chained_res$edges_df, chained_res$nodes_df)
ok = .report("chained: igraph vcount/ecount match nodes_df/edges_df",
             igraph::vcount(gr) == nrow(chained_res$nodes_df) && igraph::ecount(gr) == nrow(chained_res$edges_df)) && ok
ok = .report("chained: igraph adjacency matches TCRdist_to_sparse_matrix",
             all(as.matrix(igraph::as_adjacency_matrix(gr)) == as.matrix(adj_mat))) && ok

results["chained_adds_stay_consistent"] = ok

## ===========================================================================
## Test 4: input validation
## ===========================================================================
ok = TRUE
ok = .report("errors on edges_df missing required columns",
             inherits(tryCatch(add_to_tcr_network(data.frame(a = 1), part1_res$nodes_df, part2),
                                error = function(e) e), "error")) && ok
ok = .report("errors on nodes_df missing tcr_index",
             inherits(tryCatch(add_to_tcr_network(part1_res$edges_df, data.frame(a = 1), part2),
                                error = function(e) e), "error")) && ok
ok = .report("errors on empty nodes_new_df",
             inherits(tryCatch(add_to_tcr_network(part1_res$edges_df, part1_res$nodes_df, part2[0, ]),
                                error = function(e) e), "error")) && ok
results["input_validation"] = ok

cat("\n=== Summary ===\n")
print(results)
if (all(results)) {
  cat("\nAll tests PASSED.\n")
} else {
  stop("add_to_tcr_network() tests failed: ", paste(names(results)[!results], collapse = ", "))
}
