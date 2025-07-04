#' @title
#' Plot the number of partners for each single chain that is matched to another chain.
#'
#' @description
#' \code{plot_num_partners()} creates bar plots for alpha and beta chains showing how many
#' partners they are paired with by the MAD-HYPE and/or T-shell algorithms.
#'
#' @details
#' For each sample, the function creates stacked bar charts for alpha and beta chains,
#' showing the proportion of them (among all called pairs) that are paired with 1 chain,
#' 2 chains, 3 chains, etc.
#'
#' @param data a dataset created by \code{\link{load_tirtlseq}()} and possibly \code{\link{filter_dataset}()}
#' @param group_col a column of the metadata to use to group multiple samples into one bar plot
#' @param fraction whether to plot the fraction of chains or the total number of chains
#' (default is TRUE, i.e. plot fractions)
#' @param include_non_functional whether to include chains with non-functional
#' cdr3 sequences when tabulating the output.
#' @param max_partners the maximum number of partners, N, to include in the plots.
#' All chains with more than N partners will be grouped together under the ">N" category.
#' @param return_data if TRUE, return the data used to make the plots
#'
#' @return
#' Either a bar chart (ggplot object) with facets (sub-plots) for each sample or
#' a list with two objects:
#'
#' \code{$plot} the plot referenced above
#'
#' \code{$data} the data used to create the plot
#'
#' @seealso \code{\link{identify_non_functional_seqs}()}
#'
#' @export
#' @examples
#' # example code
#'
#'

plot_num_partners = function(data,
                             #label_col = NULL,
                             group_col = NULL,
                             fraction = TRUE,
                             include_non_functional = FALSE,
                             max_partners = 5,
                             samples = NULL,
                             return_data = FALSE,
                             color_scheme = NULL
                             ) {

  meta = data$meta
  data = data$data

  if(!is.null(samples)) {
    ind = match(samples, names(data))
    data = data[ind]
    meta = meta[ind,]
  }

  if(include_non_functional) {
    data = lapply(data, function(x) x[["paired"]]) %>% setNames(names(data))
  } else {
    data = lapply(data, function(x) x[["paired"]] %>% filter(is_functional)) %>% setNames(names(data))
  }


  is_paired = .is.paired(data)
  is_list = .is.list.only(data)
  if(!is_paired) stop("'data' must be paired chain output from TIRTL-seq")
  data = remove_dupes_paired(data)

  if(is_list) {
    gg_df = lapply(1:length(data), function(i) {
      df_tmp = .get_num_partners_single(data[[i]], max_partners = max_partners)
      df_tmp$sample_num = i
      if(!is.null(group_col)) {
        df_tmp$Group = meta[[group_col]][i]
      } else {
        df_tmp$Group = meta[[1]][i]
      }
      #if(!is.null(label_col)) df_tmp$label = meta[[label_col]][i]
      return(df_tmp)
      }) %>% bind_rows() %>% group_by(Group, n_partners, chain) %>%
      summarize(Frequency = sum(Frequency), Fraction = mean(Fraction)) ## note this is taking the mean proportion over samples rather than summing all frequencies and dividing by the total
  } else {
    gg_df = .get_num_partners_single(data)
  }
  var = ifelse(fraction, sym("Fraction"), sym("Frequency"))
  char = paste(">", max_partners, sep = "")
  lvls = c(1:max_partners, char) %>% rev()
  gg_df$n_partners = factor(gg_df$n_partners, levels = lvls)
  gg = ggplot(gg_df) +
    geom_col(aes(x = chain, y = !!var, fill = n_partners), color="black", size = 0.25) +
    xlab("") +
    theme_classic()
  gg = gg + scale_fill_manual(values = .tirtl_colors_gradient(palette=color_scheme,length(lvls) ))
  if(is_list) {
    #facet_formula = as.formula(paste("~",group_col, sep = ""))
    gg = gg + facet_wrap(~Group)
  }
  if(return_data) {
    res = list(plot=gg, data = gg_df)
  } else {
    res = gg
  }
  return(res)
}

.get_num_partners_single = function(df, max_partners = 5) {
  df = remove_dupes_paired(df)
  alpha_tbl = table(df$alpha_nuc) %>% table() %>% as.data.frame.table() %>%
    magrittr::set_colnames(c("n_partners", "Frequency")) %>%
    mutate( Fraction = Frequency/sum(Frequency), chain = "alpha" ) %>%
    mutate(n_partners = as.integer(as.character(n_partners)))
  beta_tbl = table(df$beta_nuc) %>% table() %>% as.data.frame.table() %>%
    magrittr::set_colnames(c("n_partners", "Frequency")) %>%
    mutate( Fraction = Frequency/sum(Frequency), chain = "beta" ) %>%
    mutate(n_partners = as.integer(as.character(n_partners)))
  df_long = bind_rows(alpha_tbl, beta_tbl)
  gt_char = paste(">", max_partners, sep = "")
  df_less_max = df_long %>% filter(n_partners <= max_partners) %>%
    mutate(n_partners = as.character(n_partners))
  df_greater = df_long %>% filter(n_partners > max_partners) %>%
    group_by(chain) %>%
    summarize(Frequency = sum(Frequency), Fraction = sum(Fraction)) %>%
    ungroup() %>%
    mutate(n_partners = gt_char)
  df_out = bind_rows(df_less_max, df_greater) #%>%
    #mutate(n_partners = factor(n_partners, levels = c(as.character(sort(as.integer(unique(df_less_max$n_partners)))), gt_char)) )
  return(df_out)
}
