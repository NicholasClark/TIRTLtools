#' @title
#' Calculate and plot agreement/overlap between samples
#'
#' @description
#' `sample_overlap()` returns a heatmap showing the overlap in most frequent TCRs
#' among pairs of samples in a dataset
#'
#' @details
#' The function takes the top N most frequent TCRs found in each dataset (default 200) and
#' compares their overlap between samples.
#'
#' @param data the dataset, an object loaded using the `load_tirtlseq()` function
#' @param chain which chain to plot: either paired or alpha-/beta-pseudobulk. (default "paired")
#' @param n_seq the number of most frequent TCR sequences to use (default 200)
#' @param show_row_names whether to show row names for the heatmap (default FALSE)
#' @param show_column_names whether to show column names for the heatmap (default FALSE)
#' @param label_col a column of the metadata to use as labels for rows and columns (default "Sample", uses the sample_id)
#' @param title a title for the heatmap
#'
#' @return
#' A heatmap with hierarchically clustered rows and columns showing the number of
#' TCRs shared between each pair of samples, among their top N most frequent TCRs.
#'
#' @seealso [func1()], [func2()], and [func3()] for similar functions
#'
#' @export
#' @examples
#' # example code
#'
#'

sample_overlap = function(data, chain = c("paired", "alpha", "beta"),
                          n_seq = 200,
                          show_row_names = FALSE, show_column_names = FALSE,
                          label_col = "Sample",
                          title = "") {
  meta = data$meta
  data = data$data
  chain = chain[1]

  #labels = get_labels_from_col(meta, label_col)

  df_list = lapply(data, function(x) x[[chain]]) %>% setNames(names(data))
  seq_list = lapply(df_list, function(df_tmp) {
    if("readCount" %in% colnames(df_tmp)) { ### pseudo-bulk data
      #top_seqs = df_tmp[order(readCount, decreasing = TRUE)]$targetSequences[1:n_seq]
      data.table::setkey(df_tmp, readCount)
      n_row = dim(df_tmp)[1]
      top_seqs = df_tmp$targetSequences[(n_row-n_seq+1):n_row]
    } else if("wij" %in% colnames(df_tmp)) { ### paired TIRTL-seq data
      data.table::setkey(df_tmp, wij)
      n_row = dim(df_tmp)[1]
      top_seqs = df_tmp$alpha_beta[(n_row-n_seq+1):n_row]
    }
    return(top_seqs)
  })
  olap_mat = sapply(seq_list, function(x) {
    sapply(seq_list, function(y) {
      length(intersect(x,y))
    })
  })
  if(label_col == "Sample") {
    meta = meta %>% as.data.frame() %>% tibble::column_to_rownames(colnames(meta)[1])
    meta$label = NULL
  } else {
    meta = meta %>% as.data.frame() %>% tibble::column_to_rownames(label_col)
    meta$label = NULL
  }

  labels = rownames(meta)
  rownames(olap_mat) = labels
  colnames(olap_mat) = labels

  jaccard_mat = sapply(seq_list, function(x) {
    sapply(seq_list, function(y) {
      length(intersect(x,y))/length(union(x,y))
    })
  })
  rownames(jaccard_mat) = labels
  colnames(jaccard_mat) = labels
  row_ann = ComplexHeatmap::rowAnnotation(df = meta)

  col_hc <- hclust(as.dist(1-jaccard_mat), method = "average")
  row_hc = col_hc
  #row_hc <- hclust(as.dist(1-jaccard_mat), method = "average")

  ComplexHeatmap::Heatmap(olap_mat, left_annotation = row_ann,
                          show_row_names = show_row_names,
                          show_column_names = show_column_names,
                          #clustering_distance_columns = as.dist(1-jaccard_mat),
                          #clustering_distance_rows = as.dist(1-jaccard_mat),
                          cluster_rows = dendsort::dendsort(as.dendrogram(row_hc), type="average"),
                          cluster_columns = dendsort::dendsort(as.dendrogram(col_hc), type="average"),
                          column_title = title)
}
