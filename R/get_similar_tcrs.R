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
#' @param count_by what to count when testing for enrichment of antigen annotations:
#' "vdj_db" (default) counts similar VDJdb TCRs, "submitted" counts submitted TCRs.
#' See \code{Details}.
#'
#' @details
#' For each value of "antigen.species", "antigen.gene", and "antigen.epitope" (e.g. "CMV"),
#' a one-sided test checks whether the matches between the submitted TCRs and VDJdb
#' involve that value more often than expected by chance, given how common the value is
#' among all searched VDJdb TCRs. There are two ways to count:
#'
#' * \code{count_by = "vdj_db"}: counts the similar VDJdb TCRs with the value. A
#' hypergeometric test treats the similar VDJdb TCRs as a random draw from all searched
#' VDJdb TCRs (the same test used for gene set enrichment). VDJdb contains many
#' near-identical TCRs for some epitopes, so one submitted TCR can contribute many matches.
#' * \code{count_by = "submitted"}: counts the submitted TCRs that are similar to at least
#' one VDJdb TCR with the value. Only submitted TCRs with at least one match are included.
#' If a submitted TCR matches m VDJdb TCRs, then under the null hypothesis those m TCRs are
#' a random draw from VDJdb, which gives its probability of matching the value at least
#' once. The count of submitted TCRs then follows a Poisson-binomial distribution, which
#' is used to calculate an exact p-value.
#'
#' @returns a list with four data frames:
#' * \code{similar_tcrs}: the rows of \code{\link{vdj_db}} that are similar to at least
#' one of the supplied TCRs, with two added columns: \code{n_similar} (the number of supplied
#' TCRs within the cutoff) and \code{min_TCRdist} (the smallest TCRdist to a supplied TCR).
#' * \code{antigen_species}, \code{antigen_gene}, \code{antigen_epitope}: enrichment
#' statistics for each value of the column. See \code{Details} and the columns below.
#'
#' The enrichment data frames have the columns:
#' * \code{value}: the value of the column (e.g. "CMV")
#' * \code{n_similar_db}: the number of similar VDJdb TCRs with this value
#' * \code{n_similar_db_total}: the total number of similar VDJdb TCRs
#' * \code{n_db}: the number of searched VDJdb TCRs with this value
#' * \code{n_db_total}: the total number of searched VDJdb TCRs
#' * \code{n_submitted_similar}: the number of submitted TCRs similar to at least one
#' VDJdb TCR with this value
#' * \code{n_submitted_similar_total}: the number of submitted TCRs similar to at least one
#' VDJdb TCR
#' * \code{n_submitted_total}: the total number of submitted TCRs
#' * \code{count_by}: what was counted for the test ("vdj_db" or "submitted")
#' * \code{expected}: the expected count under the null hypothesis (of \code{n_similar_db}
#' if \code{count_by = "vdj_db"}, of \code{n_submitted_similar} if \code{count_by = "submitted"})
#' * \code{fold_enrichment}: the observed count divided by \code{expected}
#' * \code{p_value}: one-sided p-value for enrichment
#' * \code{p_adj}: p-value adjusted for multiple testing (Benjamini-Hochberg)
#'
#' @family tcr_similarity
#' @seealso \code{\link{TCRdist}()}
#' @export
#'
#' @examples
#' \dontrun{
#' res = get_similar_tcrs(tcr_df)
#' res$similar_tcrs
#' res$antigen_species
#' ## count submitted TCRs instead of VDJdb TCRs
#' res = get_similar_tcrs(tcr_df, count_by = "submitted")
#' }
get_similar_tcrs = function(
    tcr_df,
    remove_MAIT = FALSE,
    params = NULL,
    submat = NULL,
    tcrdist_cutoff = NULL,
    chunk_size = 1000,
    backend = c("auto", "cupy", "mlx", "cpp"),
    count_by = c("vdj_db", "submitted")
) {
  count_by = match.arg(count_by)

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

  ## VDJ-db nodes that have an edge with at least one query TCR
  nodes_df_similar = nodes_df |>
    filter(source == "vdj_db") |>
    inner_join(edge_stats, by = c("tcr_index" = "node2_idx")) |>
    arrange(min_TCRdist)

  edge_stats_obs = edges_df |>
    group_by(node1_idx) |>
    summarize(n_similar = n(), min_TCRdist = min(TCRdist), .groups = "drop")

  ## query TCRs that have at least one hit in VDJ-db
  nodes_df_query_found = nodes_df |>
    filter(source == "observed") |>
    inner_join(edge_stats_obs, by = c("tcr_index" = "node1_idx")) |>
    arrange(min_TCRdist)

  nodes_df_return_db = nodes_df_similar[,c("n_similar", "min_TCRdist", "tcr_index", "source", colnames(db))]
  nodes_df_return_obs = nodes_df_query_found[,c("n_similar", "min_TCRdist", "tcr_index", "source", colnames(tcr_df))]

  ## enrichment of antigen annotations among the similar VDJdb TCRs
  db_nodes = nodes_df |> filter(source == "vdj_db")
  n_submitted = sum(nodes_df$source == "observed")
  stats = lapply(c("antigen.species", "antigen.gene", "antigen.epitope"), function(col) {
    .similar_tcr_enrichment(db_nodes, edges_df, column = col, n_submitted = n_submitted,
                            count_by = count_by)
  }) |> setNames(c("antigen_species", "antigen_gene", "antigen_epitope"))

  return(c(list(similar_vdj_db = nodes_df_return_db, similar_query = nodes_df_return_obs), stats))
}

#' Enrichment statistics of a VDJdb annotation column among similar TCRs
#'
#' For each value of `column` (e.g. "CMV" for "antigen.species"), a one-sided test of
#' whether matches to VDJdb TCRs with that value are more common than expected, given its
#' frequency among all searched VDJdb TCRs. Missing or empty values are not tested, but are
#' counted in the totals.
#'
#' * count_by = "vdj_db": hypergeometric test on the number of similar VDJdb TCRs.
#' * count_by = "submitted": test on the number of submitted TCRs matching the value at
#'   least once. Submitted TCR i with m_i matches hits the value with probability
#'   1 - P(no draws with the value in m_i draws from VDJdb), and the total follows a
#'   Poisson-binomial distribution.
#'
#' @param db_nodes the VDJdb nodes that were searched (with tcr_index and `column`)
#' @param edges_df oriented edges: node1_idx is a submitted TCR, node2_idx a VDJdb TCR
#' @param column the name of the annotation column in db_nodes
#' @param n_submitted the number of submitted TCRs
#' @param count_by "vdj_db" or "submitted"
#' @returns a data frame with one row per value of `column`, sorted by p-value
#' @noRd
.similar_tcr_enrichment = function(db_nodes, edges_df, column, n_submitted,
                                   count_by = c("vdj_db", "submitted")) {
  if(column == "antigen.gene") db_nodes = db_nodes |> mutate(antigen.gene = glue::glue("{antigen.species} {antigen.gene}"))
  if(column == "antigen.epitope") db_nodes = db_nodes |> mutate(antigen.epitope = glue::glue("{antigen.species} {antigen.gene} {antigen.epitope}"))

  count_by = match.arg(count_by)
  values = as.character(db_nodes[[column]])
  values[values == ""] = NA
  edges_df = edges_df |> distinct(node1_idx, node2_idx)
  is_similar = db_nodes$tcr_index %in% edges_df$node2_idx

  N = nrow(db_nodes) ## all searched VDJdb TCRs
  n = sum(is_similar) ## similar VDJdb TCRs

  ## number of submitted TCRs with at least one similar VDJdb TCR of each value
  edge_values = values[match(edges_df$node2_idx, db_nodes$tcr_index)]
  submitted_counts = tapply(edges_df$node1_idx, edge_values, function(x) length(unique(x)))
  ## number of VDJdb matches for each submitted TCR with at least one match
  n_matches = as.vector(table(edges_df$node1_idx))

  res = tibble(value = values, is_similar = is_similar) |>
    filter(!is.na(value)) |>
    group_by(value) |>
    summarize(n_similar_db = sum(is_similar), n_db = n(), .groups = "drop") |>
    mutate(
      n_similar_db_total = n,
      n_db_total = N,
      n_submitted_similar = as.integer(submitted_counts[value]),
      n_submitted_similar = ifelse(is.na(n_submitted_similar), 0L, n_submitted_similar),
      n_submitted_similar_total = length(n_matches),
      n_submitted_total = n_submitted,
      count_by = count_by
    )

  if(count_by == "vdj_db") {
    res = res |> mutate(
      expected = n * n_db / N,
      fold_enrichment = n_similar_db / expected,
      p_value = stats::phyper(n_similar_db - 1, n_db, N - n_db, n, lower.tail = FALSE)
    )
  } else {
    stats_mat = t(vapply(seq_len(nrow(res)), function(i) {
      ## chance that each matched submitted TCR hits this value at least once
      p_hit = 1 - stats::phyper(0, res$n_db[i], N - res$n_db[i], n_matches)
      k = res$n_submitted_similar[i]
      c(expected = sum(p_hit), p_value = if(k == 0) 1 else .poisson_binomial_upper(k, p_hit))
    }, numeric(2)))
    res = res |> mutate(
      expected = stats_mat[, "expected"],
      fold_enrichment = n_submitted_similar / expected,
      p_value = stats_mat[, "p_value"]
    )
  }

  res |>
    mutate(p_adj = stats::p.adjust(p_value, method = "BH")) |>
    select(value, n_similar_db, n_similar_db_total, n_db, n_db_total,
           n_submitted_similar, n_submitted_similar_total, n_submitted_total,
           count_by, expected, fold_enrichment, p_value, p_adj) |>
    arrange(p_value, desc(fold_enrichment))
}

#' Upper tail P(X >= k) of a Poisson-binomial distribution
#'
#' X is the number of successes in independent trials with success probabilities `p`,
#' computed exactly by dynamic programming over the trials.
#' @noRd
.poisson_binomial_upper = function(k, p) {
  if(k > length(p)) return(0)
  dens = 1 ## dens[j + 1] = P(X = j) after the trials so far
  for(p_i in p) dens = c(dens * (1 - p_i), 0) + c(0, dens * p_i)
  min(1, sum(dens[(k + 1):length(dens)]))
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
