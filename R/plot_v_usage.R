#' @title
#' Plotting of V-segment usage
#'
#' @description
#' `plot_v_usage` returns a bar plot of the usage of VA or VB segments.
#'
#' @details
#' This function ...
#'
#' @param df a data frame containing paired TCRs from TIRTL-seq output
#' @param chain the chain of the TCRs to plots - "alpha" for VA segments or "beta" for VB segments
#' @param n_max the maximum number of V-segments to plot. The top `n_max` most common V-segments will be plotted
#'
#' @returns a ggplot object with a bar chart of the most common V-segments.
#'
#' @examples
#'

plot_v_usage = function(df, chain = c("alpha", "beta"), n_max = 25) {
  chain = match.arg(chain)
  if(!chain %in% c("alpha", "beta")) stop()
  if(chain == "alpha") var = "va"
  if(chain == "beta") var = "vb"
  df_summ = df %>% group_by(!!sym(var)) %>% summarize(n = n()) %>% arrange(desc(n))
  #df_gg = df_summ[1:n_max,] %>% mutate(va = factor(!!var, levels = rev(!!var)))
  df_gg = df_summ[1:n_max,]
  df_gg[[var]] = factor(df_gg[[var]], levels = rev(df_gg[[var]]))

  xlab_title = "Count"
  ggplot(df_gg) + geom_col(aes(x=n, y=!!sym(var), fill = n), color = "grey20") + theme_classic() +
    ylab("") + xlab(xlab_title) +
    paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", direction = -1, name = xlab_title)

  # df_summ = df %>% group_by(va) %>% summarize(n = n()) %>% arrange(desc(n))
  # df_gg = df_summ[1:n_max,] %>% mutate(va = factor(va, levels = rev(va)))
  # ggplot(df_gg) + geom_col(aes(x=n, y=va)) + theme_classic()
}
