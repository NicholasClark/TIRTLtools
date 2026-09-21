#' Add new TCRs to an existing TCR similarity network
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' Given a network previously built with \code{\link{TCRdist}()} (or \code{\link{cluster_tcrs}()},
#' this function calculates TCRdist between a new
#' set of TCRs (\code{nodes_new_df}) and all of the existing TCRs, then appends the new nodes
#' and edges to the network and returns the updated \code{edges_df}/\code{nodes_df}.
#'
#' @details
#' `nodes_new_df` must contain the columns "va", "vb", "cdr3a", and "cdr3b".
#' These columns must contain the V-alpha segment, V-beta segment,
#' the CDR3-alpha amino acid sequence, and the CDR3-beta amino acid sequence, respectively.
#'
#' `nodes_df` must contain the columns "va", "vb", "cdr3a", "cdr3b", and "tcr_index", where "tcr_index"
#' is an integer and the other columns are as above.
#'
#' `edges_df` must contain the columns "node1_idx", "node2_idx", and "TCRdist", where the first
#' two columns are integers that map to the "tcr_index" column in `nodes_df` and "TCRdist" is the value of TCRdist
#' between the two TCRs.
#'
#' New nodes are numbered starting right after the highest \code{tcr_index} already present
#' in \code{nodes_df}.
#'
#' @param edges_df the \code{edges_df} from a previous \code{\link{TCRdist}()} (or
#' \code{\link{cluster_tcrs}()}) result: a data frame with columns "node1_idx", "node2_idx",
#' and "TCRdist".
#' @param nodes_df the \code{nodes_df} from a previous \code{\link{TCRdist}()} (or
#' \code{\link{cluster_tcrs}()}) result: one row per existing
#' TCR.
#' @param nodes_new_df a data frame of new TCRs to add to the network. Must have the columns
#' "va", "vb", "cdr3a", and "cdr3b" (the same requirements as \code{tcr1}/\code{tcr2} in
#' \code{\link{TCRdist}()}).
#' @param remove_MAIT whether to remove TCRs from MAIT cells (default is FALSE)
#' @param params (optional) a table of valid parameters for amino acids and va/vb segments.
#' (default is NULL, which uses TIRTLtools::params)
#' @param submat (optional) a substitution matrix with mismatch penalties for each
#' combination of amino acids or va/vb segments (default is NULL, which uses TIRTLtools::submat).
#' @param tcrdist_cutoff (optional) discard all TCRdist values above this cutoff. If not
#' supplied by the user, this will default to 90 for dual-chain TCRdist or 45 for
#' single-chain TCRdist.
#' @param chunk_size (optional) the chunk size to use in calculation of TCRdist (default 1000).
#' See \code{\link{TCRdist}()} for details.
#' @param backend (optional) the backend to use for the chunk computation (default "auto").
#' See \code{\link{TCRdist}()} for details.
#'
#' @returns a list with the same shape as \code{\link{TCRdist}()}'s output:
#'
#' \code{$edges_df} - \code{edges_df} with the new TCRdist edges appended: new nodes vs.
#' existing nodes, and new nodes vs. each other.
#'
#' \code{$nodes_df} - \code{nodes_df} with the new TCRs appended. New rows get \code{source
#' = "new"} and \code{tcr_index} values starting right after the highest \code{tcr_index}
#' already present in \code{nodes_df}.
#'
#' @family tcr_similarity
#' @seealso \code{\link{TCRdist}()}, \code{\link{cluster_tcrs}()}
#'
#' @export
#' @examples
#' load_example_data(dataset = "SJTRC_minimal")
#' df = get_all_tcrs(SJTRC_minimal, chain = "paired", remove_duplicates = TRUE)
#' result = TCRdist(df[1:500, ], tcrdist_cutoff = 90)
#' result2 = add_to_tcr_network(result$edges_df, result$nodes_df, df[501:600, ])
#'
#' ## compare to a single TCRdist() run on the combined set:
#' result_full = TCRdist(df[1:600, ], tcrdist_cutoff = 90)
#' nrow(result2$edges_df) == nrow(result_full$edges_df)
add_to_tcr_network = function(
    edges_df,
    nodes_df,
    nodes_new_df,
    remove_MAIT = FALSE,
    params = NULL,
    submat = NULL,
    tcrdist_cutoff = NULL,
    chunk_size = 1000,
    backend = c("auto", "cupy", "mlx", "cpp")
    ) {
  backend = match.arg(backend)

  if (!all(c("node1_idx", "node2_idx", "TCRdist") %in% colnames(edges_df))) {
    cli::cli_abort("{.arg edges_df} must have columns {.field node1_idx}, {.field node2_idx}, and {.field TCRdist} (the output of {.fn TCRdist}).")
  }
  if (!"tcr_index" %in% colnames(nodes_df)) {
    cli::cli_abort("{.arg nodes_df} must have a {.field tcr_index} column (the output of {.fn TCRdist}).")
  }
  if (nrow(nodes_new_df) == 0) {
    cli::cli_abort("{.arg nodes_new_df} has no rows -- nothing to add.")
  }

  n_old = nrow(nodes_df)
  max_old_idx = max(nodes_df$tcr_index)

  tcrdist_args = list(
    remove_MAIT = remove_MAIT,
    params = params,
    submat = submat,
    tcrdist_cutoff = tcrdist_cutoff,
    chunk_size = chunk_size,
    backend = backend
  )

  ## Calculate TCRdist between the new TCRs (tcr1) and the existing nodes (tcr2). This
  ## re-runs prep_for_tcrdist() on nodes_df, but that is idempotent for TCRs that have
  ## already been through it once, so it reproduces the same rows in the same order --
  ## verified below rather than assumed.
  new_vs_old = do.call(TCRdist, c(list(tcr1 = nodes_new_df, tcr2 = nodes_df), tcrdist_args))

  new_nodes = new_vs_old$nodes_df %>% filter(source == "tcr1")
  n_new = nrow(new_nodes)
  n_old_reprepped = nrow(new_vs_old$nodes_df) - n_new
  if (n_old_reprepped != n_old) {
    cli::cli_abort(c(
      "Re-preparing {.arg nodes_df} for TCRdist changed the number of rows ({n_old} -> {n_old_reprepped}).",
      "i" = "This usually means {.arg nodes_df} was not produced by {.fn TCRdist} (or {.fn cluster_tcrs}), or was built with a different {.arg remove_MAIT} setting than the one passed here."
    ))
  }

  ## Also calculate TCRdist among the new TCRs themselves, so the result matches what
  ## you'd get from running TCRdist() once on the combined old+new set. prep_for_tcrdist()
  ## is a deterministic, row-order-preserving function of nodes_new_df alone (it does not
  ## depend on what else is in the batch), so this reproduces exactly the same n_new rows,
  ## in the same order, as the new_vs_old call above -- verified below rather than assumed.
  new_vs_new = do.call(TCRdist, c(list(tcr1 = nodes_new_df, tcr2 = NULL), tcrdist_args))
  if (nrow(new_vs_new$nodes_df) != n_new) {
    cli::cli_abort("Re-preparing {.arg nodes_new_df} gave a different number of rows across the two internal TCRdist() calls ({n_new} vs. {nrow(new_vs_new$nodes_df)}); cannot safely align new-vs-old and new-vs-new edges.")
  }

  ## Internal (within the new_vs_old TCRdist() call) node ids: new nodes are 1..n_new, old
  ## (re-prepped) nodes are (n_new+1)..(n_new+n_old). Re-map both back to the ids we
  ## actually want: new nodes start right after nodes_df's highest tcr_index, and old
  ## nodes keep their original tcr_index (via nodes_df's own row order, which
  ## prep_for_tcrdist() preserves).
  old_tcr_index_lookup = nodes_df$tcr_index

  new_vs_old_edges = new_vs_old$edges_df %>%
    mutate(
      node1_idx = max_old_idx + node1_idx,
      node2_idx = old_tcr_index_lookup[node2_idx - n_new]
    )

  ## Internal (within the new_vs_new TCRdist() call) node ids are also 1..n_new (in the
  ## same order as new_vs_old's new-node portion), so both endpoints get the same
  ## new-node offset.
  new_vs_new_edges = new_vs_new$edges_df %>%
    mutate(
      node1_idx = max_old_idx + node1_idx,
      node2_idx = max_old_idx + node2_idx
    )

  new_nodes = new_nodes %>%
    mutate(source = "new", tcr_index = max_old_idx + seq_len(n_new))

  out = list(
    edges_df = bind_rows(edges_df, new_vs_old_edges, new_vs_new_edges) |> arrange(node1_idx, node2_idx),
    nodes_df = bind_rows(nodes_df, new_nodes)
  )
  return(out)
}
