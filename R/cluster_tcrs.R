#' @title
#' Clustering of TCRs using TCRdist and the Leiden algorithm
#'
#' @description
#' The \code{cluster_tcrs()} function aggregates all of the paired TCRs found in the data,
#' calculates pairwise similarity using the va, vb, cdr3a, and cdr3b regions (via TCRdist),
#' and clusters the results using the Leiden algorithm.
#'
#' @details
#' The function also filters the dataset to TCRs that are valid for TCRdist.
#'
#' The following TCRs are removed:
#'  - TCRs that contain stop codons (*) or frame shifts (_) in their cdr3a or cdr3b regions
#'  - TCRs that contain a cdr3 region with 5 or less amino acids
#'  - TCRs that contain a v segment allele not found in our parameter table
#'
#' V-segments that do not specify an allele (e.g. "TRAV1-2" instead of "TRAV1-2*01")
#' will be assigned to the "*01" allele.
#'
#' @param data a list of TIRTLseq TCR data for samples created with \code{\link{load_tirtlseq}()}
#' @param tcrdist_cutoff the \code{\link{TCRdist}()} function will only record TCRdist values
#' less than or equal to the cutoff. Default is 90. Note: Higher cutoffs will return more
#' data, at most NxN where N is the number of unique TCRs.
#' @param resolution the "resolution" parameter for the Leiden algorithm. A lower
#' value will produce larger clusters and a higher value will produce smaller clusters.
#' Typical values are in the 0.1 - 2.0 range. A higher value may be better for densely
#' connected data while a lower value may be better for moderately connected data. Default is 0.1.
#' @param with_vdjdb if TRUE, observed clones will be compared and clustered with annotated clones from VDJ-db.
#' If parameter is a data frame, the supplied data frame will be used as the database.
#' @param allow_self_edges (default TRUE) if FALSE, only calculate TCRdist between members of the input data and vdj_db
#'
#' @return
#' Returns a list with the following elements:
#'
#' \code{$df} - a data frame with all unique TCRs along with cluster annotations
#'
#' \code{$dist_df} - a data frame with distances (TCRdist) between TCR pairs in long format
#'
#' \code{$sparse_adj_mat} - an adjacency matrix (in sparse format) marking TCR pairs with TCRdist <= tcrdist_cutoff
#'
#' \code{$graph_adj} - an igraph object created from the adjacency matrix
#'
#' \code{$tcrdist_cutoff} - the cutoff used for TCRdist
#'
#' \code{$resolution} - the resolution parameter used for the Leiden algorithm
#'
#' @seealso \code{\link{plot_clusters}()}, \code{\link{identify_non_functional_seqs}()}, \code{\link{TCRdist}()}
#'
#' @export
#' @examples
#' # example code
#' # paired = load_tirtlseq("your_directory/")
#' # obj = cluster_tcrs(paired)
#'
cluster_tcrs = function(data, tcrdist_cutoff = 90, resolution = 0.1, with_vdjdb = TRUE, allow_self_edges = TRUE) {
  cluster.type = "leiden"
  chain = "paired"
  if(.is.list.only(data)) { ### if an object loaded by load_tirtlseq
    df_all_obs = get_all_tcrs(data, chain, remove_duplicates = TRUE) %>%
      mutate(source = "observed") ## get all tcrs in one data frame
  } else { ### if a dataframe with TCRs
    df_all_obs = data %>% mutate(source = "observed")
  }
  if(is.data.frame(with_vdjdb)) {
    vdj = with_vdjdb
    df_all = bind_rows(df_all_obs, vdj)
  } else {
    if(with_vdjdb) {
      vdj = TIRTLtools::vdj_db %>% mutate(source = "vdj-db")
      df_all = bind_rows(df_all_obs, vdj)
    } else {
      df_all = df_all_obs
    }
  }

  if( (!allow_self_edges) && with_vdjdb ) {
    dist = TCRdist(tcr1 = df_all_obs, tcr2 = vdj, tcrdist_cutoff = tcrdist_cutoff)
  } else {
    dist = TCRdist(df_all, tcrdist_cutoff = tcrdist_cutoff)
  }
  dist_df = dist$TCRdist_df
  dist_input = dist$tcr1

  dist_df = dist_df %>% mutate(
    edge1_1index = edge1_0index + 1,
    edge2_1index = edge2_0index + 1
  )
  dist_df$weight_binary = 1
  #dist_df$TCRdist_mod = ifelse(dist_df$TCRdist_mod == 0, -1, dist_df$TCRdist_mod)
  n_valid = dim(dist_input)[1]
  sparse_weight_mat_binary = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
                                                  x=dist_df$weight_binary, symmetric = TRUE,
                                                  dims = c(n_valid, n_valid))
  # sparse_tcrdist_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
  #                                           x=dist_df$TCRdist_mod, symmetric = TRUE,
  #                                           dims = c(n_valid, n_valid))

  # sparse_weight_mat = Matrix::sparseMatrix(i=dist_df$edge1_1index, j=dist_df$edge2_1index,
  #                                          x=dist_df$weight, symmetric = TRUE,
  #                                          dims = c(n_valid, n_valid) )

  gr_binary = igraph::graph_from_adjacency_matrix(sparse_weight_mat_binary, mode = "undirected", weighted = NULL)
  leiden_clust = igraph::cluster_leiden(gr_binary, resolution = resolution)
  dist_input$idx_1index = 1:dim(dist_input)[1]
  dist_input$cluster = leiden_clust$membership ### assign clusters to TCR sequences
  tab = table(leiden_clust$membership)
  single = sum(tab == 1)
  g2 = sum(tab >= 2)
  g10 = sum(tab >= 10)
  g50 = sum(tab >= 50)
  g100 = sum(tab >= 100)
  msg1 = paste("Out of ", dim(dist_input)[1], " valid TCRs, ", g2, " clusters detected and ", single, " singleton TCRs.", sep = "") %>% .add_newline()
  msg2 = paste(g10, " clusters of size >= 10, ", g50, " clusters of size >= 50, ", g100, " clusters of size >=100.", sep = "") %>% .add_newline()
  cat(msg1); cat(msg2)
  out = list(
    df = dist_input,
    dist_df = dist_df,
    #sparse_tcrdist_mat = sparse_tcrdist_mat, ### returning this as a sparse matrix is problematic because missing entries will be seen as TCRdist = 0 instead of TCRdist > cutoff.
    sparse_adj_mat = sparse_weight_mat_binary,
    graph_adj = gr_binary,
    tcrdist_cutoff = tcrdist_cutoff,
    resolution = resolution
    )
  return(out)
}
