#' Find TCRs in VDJdb that are similar to a set of TCRs
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' This function calculates TCRdist between the supplied TCRs and the TCRs in
#' \code{\link{vdj_db}} and returns the VDJdb TCRs that are within the TCRdist cutoff of
#' at least one supplied TCR.
#'
#' @param tcr_df a data frame of TCRs, in the same format as the input to \code{\link{TCRdist}()}
#' @param remove_MAIT whether to remove MAIT cells before calculating TCRdist (default is FALSE)
#' @param params (optional) a table of valid parameters for amino acids and va/vb segments
#' (default is NULL, which uses \code{\link{params}})
#' @param submat (optional) a substitution matrix with mismatch penalties
#' (default is NULL, which uses \code{\link{submat}})
#' @param tcrdist_cutoff (optional) the maximum TCRdist for two TCRs to be considered
#' similar (default is NULL, which uses 90 for paired TCRs and 45 for single chains)
#' @param chunk_size the number of rows to compare at once (default is 1000)
#' @param backend the backend to use for TCRdist calculations: "auto", "cupy", "mlx", or "cpp"
#' (default is "auto"). See \code{\link{TCRdist}()} for details.
#'
#' @returns a data frame with the rows of \code{\link{vdj_db}} that are similar to at least
#' one of the supplied TCRs, with two added columns: \code{n_similar} (the number of supplied
#' TCRs within the cutoff) and \code{min_TCRdist} (the smallest TCRdist to a supplied TCR).
#'
#' @family tcr_similarity
#' @seealso \code{\link{TCRdist}()}
#' @export
#'
#' @examples
#' \dontrun{
#' similar = get_similar_tcrs(tcr_df)
#' }
get_similar_tcrs = function(
    tcr_df,
    remove_MAIT = FALSE,
    params = NULL,
    submat = NULL,
    tcrdist_cutoff = NULL,
    chunk_size = 1000,
    backend = c("auto", "cupy", "mlx", "cpp")
) {

  tcr_df = tcr_df |> add_ids_to_paired_df()
  tcr_df = tcr_df |> filter(!duplicated(receptor))

  db = TIRTLtools::vdj_db #|> add_ids_to_paired_df()
  td = forward_args(TCRdist, tcr1 = tcr_df, tcr2 = db)
  nodes_df = td$nodes_df
  nodes_df$source[nodes_df$source == "tcr1"] = "observed"
  nodes_df$source[nodes_df$source == "tcr2"] = "vdj_db"

  ## keep only observed <-> vdj_db edges, with node1 always observed and node2 always vdj_db
  edges_df = .orient_edges(td$edges_df, nodes_df, from = "observed", to = "vdj_db")

  edge_stats = edges_df |>
    group_by(node2_idx) |>
    summarize(n_similar = n(), min_TCRdist = min(TCRdist), .groups = "drop")

  nodes_df_similar = nodes_df |>
    filter(source == "vdj_db") |>
    inner_join(edge_stats, by = c("tcr_index" = "node2_idx")) |>
    select(-tcr_index, -source) |>
    arrange(min_TCRdist)

  nodes_df_return = nodes_df_similar[,c("n_similar", "min_TCRdist", colnames(db))]

  return(nodes_df_return)
}

#' Orient edges so that node1 is always from one source and node2 from another
#'
#' Edges between nodes of the same source, or involving any other source, are dropped.
#' Edges going the "wrong" way (node1 from `to`, node2 from `from`) are flipped by
#' swapping every pair of "node1_*"/"node2_*" columns (e.g. node1_idx/node2_idx, and any
#' node metadata columns added by .add_node_metadata_to_edges()).
#'
#' @param edges_df an edges data frame with node1_idx and node2_idx columns
#' @param nodes_df a nodes data frame with tcr_index and source columns
#' @param from the source that node1 should come from
#' @param to the source that node2 should come from
#' @returns the filtered, oriented edges data frame
#' @noRd
.orient_edges = function(edges_df, nodes_df, from = "observed", to = "vdj_db") {
  source1 = nodes_df$source[match(edges_df$node1_idx, nodes_df$tcr_index)]
  source2 = nodes_df$source[match(edges_df$node2_idx, nodes_df$tcr_index)]

  keep = (source1 %in% from & source2 %in% to) | (source1 %in% to & source2 %in% from)
  flip = (source1 %in% to & source2 %in% from)[keep]
  edges_df = edges_df[keep, , drop = FALSE]

  cols1 = grep("^node1_", colnames(edges_df), value = TRUE)
  cols2 = sub("^node1_", "node2_", cols1)
  paired = cols2 %in% colnames(edges_df)
  for(k in which(paired)) {
    v1 = edges_df[[cols1[k]]]
    v2 = edges_df[[cols2[k]]]
    edges_df[[cols1[k]]][flip] = v2[flip]
    edges_df[[cols2[k]]][flip] = v1[flip]
  }
  edges_df
}
