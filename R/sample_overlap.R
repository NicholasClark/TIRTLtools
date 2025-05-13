sample_overlap = function(df_list, meta, n_seq = 200,
                          show_row_names = FALSE, show_column_names = FALSE,
                          label_col = "filename",
                          title = "") {
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
  meta = meta %>% as.data.frame() %>% tibble::column_to_rownames(label_col)
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
