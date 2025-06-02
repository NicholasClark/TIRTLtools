#' @title
#' Plotting of V- or J- gene usage
#'
#' @description
#' \code{plot_gene_usage()} creates a bar chart of gene segment usage for the most common
#' V- or J- genes.
#'
#' @details
#' Fill in later, after updating/finalizing function.
#'
#' @param data either a single data frame containing paired TCRs from TIRTL-seq output or a list of
#' data frames for many experiments
#' @param groups (optional) when data is a list of data frames, this may be a vector of the same
#' length as data with labels to group the experiments.
#' @param gene the gene of the TCRs to plots
#' @param n_max the maximum number of V-segments to plot. The top \code{n_max} most common V-segments will be plotted
#'
#' @return
#' a ggplot object with a bar chart of the most common V- or J- segments.
#'
#' @export
#' @examples
#' # example here

plot_gene_usage = function(data, chain = c("paired", "alpha", "beta"),
                           group_col=NULL, gene = c("va", "vb", "ja", "jb", "v","j"), n_max = 25,
#                           groups = NULL, gene = c("va", "vb", "ja", "jb", "v","j"), n_max = 25,
                        value_type = c("auto","readFraction", "readCount", "n_wells", "readCount_max", "readCount_median", "avg", "n"),
                        agg_func = c("sum", "max"),
                        log_scale = FALSE,
                        style = c("barplot", "boxplot")) {
  meta = data$meta
  data = data$data
  chain = chain[1]
  if(!is.null(group_col)) {
    groups = meta[[group_col]]
  } else {
    groups = NULL
  }

  #labels = .get_labels_from_col(meta, label_col)

  data = lapply(data, function(x) x[[chain]]) %>% setNames(names(data))

  gene = gene[1]
  #gene = match.arg(gene)
  style = match.arg(style)
  value_type = match.arg(value_type)
  agg_func = get(agg_func[1])

  is_data_frame = is.data.frame(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])

  if(value_type == "auto") {
    if(is_paired) {
      value_type = "n_wells"
    } else {
      value_type = "readFraction"
    }
  }

  ### get top genes
  if(is_list) {
    data = lapply(1:length(data), function(i) {
      data_tmp = data[[i]]
      if(!is.null(groups)) {
        data_tmp$group = groups[i]
      } else {
        data_tmp$group = as.character(i)
      }
      return(data_tmp)
    }) %>% dplyr::bind_rows()
  }
  if(is_paired) data$n_wells = data$wij

  top_gene_df = data %>% group_by(!!sym(gene)) %>% summarize(n = n(), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE)) %>% arrange(desc(!!sym(value_type)))
  if(n_max<=dim(top_gene_df)[1]) {
    top_genes = top_gene_df[[gene]][1:n_max]
  } else {
    top_genes = top_gene_df[[gene]]
  }
  data[[gene]] = ifelse(data[[gene]] %in% top_genes, data[[gene]], "other")
  ### summarize by gene
  if(is_list && !is.null(groups)) {
    data_summ = data %>% group_by(!!sym(gene), group) %>% summarize(n = n(), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE)) %>% arrange(desc(!!sym(value_type)))
  }
  else {
    data_summ = data %>% group_by(!!sym(gene)) %>% summarize(n = n(), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE)) %>% arrange(desc(!!sym(value_type)))
  }
  data_summ[[gene]] = factor(data_summ[[gene]], levels = rev(c(top_genes, "other")))
  ### plot
  if(is_list && !is.null(groups)) {
    gg = ggplot(data_summ) + geom_col(aes(x=!!sym(value_type), y = !!sym(gene), fill = group), position = "dodge", color = "grey20")
  } else {
    gg = ggplot(data_summ) + geom_col(aes(x=!!sym(value_type), y = !!sym(gene), fill = !!sym(value_type)), color = "grey20") +
      paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1)
  }
  gg = gg + ylab("") +theme_classic()
  if(log_scale) {
    gg = gg + scale_x_log10(labels = scales::label_log())
  } else {
    gg = gg + scale_x_continuous(labels = scales::label_number_auto())
  }
  return(gg)

  # ### for list of data frames input -----------------
  # if(is_list) {
  #
  #   data_all = lapply(1:length(data), function(i) {
  #     data_tmp = data[[i]]
  #     if(!is.null(group)) data_tmp$group = groups[i]
  #     return(data_tmp)
  #   }) %>% dplyr::bind_rows()
  #
  #   if(is_paired) {
  #     ### choose v-segment to keep: top n_max most common
  #     data_summ_all = data_all %>% group_by(!!sym(gene)) %>% summarize(n = n()) %>% arrange(desc(n))
  #     v_kept = data_summ_all[[gene]][1:n_max]
  #
  #     data_all_sub = data_all %>% filter(!!sym(gene) %in% v_kept)
  #     data_gg = data_all_sub %>% group_by(group,!!sym(gene)) %>% summarize(n = n())
  #     data_gg[[gene]] = factor(data_gg[[gene]], levels = rev(v_kept))
  #
  #     #data_all_sub[[gene]] = factor(data_all_sub[[gene]], levels = rev(v_kept))
  #
  #     xlab_title = "Count"
  #     if(style == "barplot") {
  #       gg = ggplot(data_gg) + geom_col(aes(x=n, y=va, fill = group),
  #                                       color = "grey20", position = "dodge") + theme_classic() +
  #         ylab("") + xlab(xlab_title)
  #     }
  #     if(style == "boxplot") {
  #       gg = ggplot(data_gg, aes(x=n, y=va)) + geom_boxplot() +
  #         geom_point(aes(color = group)) + theme_classic() +
  #         ylab("") + xlab(xlab_title)
  #     }
  #   } else { ### pseudobulk
  #
  #   }
  # }
  # ### for single data frame input ----------------------
  # if(is_data_frame) {
  #   if(is_paired) {
  #     data$n_wells = data$wij
  #     ### sum over multiple rows w/ the same gene
  #     data_summ = data %>% group_by(!!sym(gene)) %>% summarize(n = n(), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE)) %>% arrange(desc(!!sym(value_type)))
  #
  #
  #     data_gg = data_summ[1:n_max,]
  #     data_gg[[gene]] = factor(data_gg[[gene]], levels = rev(data_gg[[gene]]))
  #
  #     xlab_title = "Count"
  #     gg = ggplot(data_gg) + geom_col(aes(x=n, y=!!sym(gene), fill = n), color = "grey20") + theme_classic() +
  #       ylab("") +
  #       paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1)
  #   } else { ### pseudobulk
  #     ### sum over multiple rows w/ the same gene
  #     data_summ = data %>% group_by(!!sym(gene)) %>% summarize(n = n(), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE)) %>% arrange(desc(!!sym(value_type)))
  #     n_genes_obs = dim(data_summ)[1]
  #     if(n_max >= n_genes_obs) {
  #       n_max = n_genes_obs
  #       keep_genes = data_summ[[gene]][1:n_max]
  #       data_gg = data_summ
  #       data_gg[[gene]] = factor(data_gg[[gene]], levels = rev(keep_genes))
  #       warning(paste('n_max >= number of genes observed: displaying', n_genes_obs, 'genes'))
  #     } else {
  #       keep_genes = data_summ[[gene]][1:n_max]
  #       data_other = data_summ[(n_max+1):dim(data_summ)[1],] %>%
  #         summarize(!!sym(gene) := "other", n=sum(n,na.rm=TRUE), !!sym(value_type) := agg_func(!!sym(value_type), na.rm = TRUE))
  #       data_gg = bind_rows(data_summ[1:n_max,], data_other)
  #       data_gg[[gene]] = factor(data_gg[[gene]], levels = c("other",rev(keep_genes)))
  #     }
  #     gg = ggplot(data_gg) + geom_col(aes(x=!!sym(value_type), y = !!sym(gene), fill = !!sym(value_type)), color = "grey20") +
  #       ylab("") +
  #       theme_classic() +
  #       paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1)
  #     if(log_scale) {
  #       gg = gg + scale_x_log10(labels = scales::label_log())
  #     } else {
  #       gg = gg + scale_x_continuous(labels = scales::label_number_auto())
  #     }
  #   }
  # }
  #return(gg)
}


# plot_v_usage_list = function(data_list, meta, gene = c("alpha", "beta"),
#                              meta_col = "patient", n_max = 25) {
#   gene = match.arg(gene)
#   if(is.data.frame(data_list)) {
#     stop("Input needs to be a list of data frames")
#   }
#   if(!gene %in% c("alpha", "beta")) stop()
#   if(gene == "alpha") gene = "va"
#   if(gene == "beta") gene = "vb"
#
#   data_all = lapply(1:length(data_list), function(i) {
#     data_tmp = data_list[[i]]
#     data_tmp$index = i
#     return(data_tmp)
#   }) %>% dplyr::bind_rows()
#   meta_cols = colnames(meta)
#   meta$index = 1:dim(meta)[1]
#   data_all = data_all %>% left_join(meta, by = "index")
#
#   ### choose v-segment to keep: top n_max most common
#   data_summ_all = data_all %>% group_by(!!sym(gene)) %>% summarize(n = n()) %>% arrange(desc(n))
#   #data_gg = data_summ[1:n_max,] %>% mutate(va = factor(!!gene, levels = rev(!!gene)))
#   v_kept = data_summ_all[[gene]][1:n_max]
#   #data_gg = data_summ_all[1:n_max,]
#   #data_gg[[gene]] = factor(data_gg[[gene]], levels = rev(data_gg[[gene]]))
#
#   data_all_sub = data_all %>% filter(!!sym(gene) %in% v_kept)
#   data_all_sub[[gene]] = factor(data_all_sub[[gene]], levels = rev(v_kept))
#
#   xlab_title = "Count"
#   ggplot(data_all_sub) + geom_bar(aes(y=!!sym(gene), fill = !!sym(meta_col)),
#                                 color = "grey20", position = "dodge") + theme_classic()# +
#     ylab("") + #xlab(xlab_title) +
#     paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1, name = xlab_title)
#
#   # data_summ = data %>% group_by(va) %>% summarize(n = n()) %>% arrange(desc(n))
#   # data_gg = data_summ[1:n_max,] %>% mutate(va = factor(va, levels = rev(va)))
#   # ggplot(data_gg) + geom_col(aes(x=n, y=va)) + theme_classic()
# }
