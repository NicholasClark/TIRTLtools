#' Make a step plot showing the extent of multi-pairing in αβ-TCRs
#'
#' @description
#' `r lifecycle::badge('experimental')`
#'
#' For this function, each αβ-TCR pair in the paired data frame is categorized by the
#' number of partners that one of the chains has (default "beta" chain). Then the fraction
#' of total pairs is plotted (y-axis) vs. the maximum number of partners allowed (x-axis).
#'
#' If `cumulative = FALSE` the fraction of pairs is plotted for each "number of partners"
#' (i.e. a probability mass function rather than a cumulative mass function).
#'
#' If `freq = TRUE`, the number of pairs rather than the fraction is plotted.
#'
#' @param data The dataset to use for plotting. Either a "TIRTLseqDataSet" containing
#' multiple samples or a "TIRTLseqData" object for one sample.
#' @param chain The chain to plot: "alpha" or "beta" or "both" (default is "beta")
#' @param cumulative Whether to plot the cumulative number of pairs (default is TRUE)
#' @param freq Whether to plot the frequency, as opposed to the fraction (default is FALSE)
#' @param step The direction of stairs (to pass to `geom_step`):  'vh' for vertical then
#' horizontal, 'hv' for horizontal then vertical.
#' @param max_x The maximum x-axis value (number of partners) allowable on the plot (default is Inf).
#'
#' @return
#' A step plot (ggplot object) as defined above.
#'
#' @family qc
#'
#' @examples
#' load_example_data(dataset = "SJTRC_longitudinal")
#'
#' plot_n_paired_step(SJTRC_longitudinal)
#' plot_n_paired_step(SJTRC_longitudinal, cumulative = FALSE)
#'
#' plot_n_paired_step(SJTRC_longitudinal, freq = TRUE)
#'
#' plot_n_paired_step(SJTRC_longitudinal, chain = "alpha")
#' plot_n_paired_step(SJTRC_longitudinal, chain = "both")
#'

plot_n_paired_step = function(
    data,
    chain = c("beta", "alpha", "both"),
    cumulative = TRUE,
    freq = FALSE,
    step = c("vh", "hv"),
    max_x = Inf
    ) {
  chain = rlang::arg_match(chain)
  step = rlang::arg_match(step)
  assert_numeric(max_x)
  assert_choice(chain, c("both","beta", "alpha"))
  assert_choice(step, c("vh", "hv"))
  assert_flag(cumulative)
  assert_flag(freq)
  assert_multi_class(data, classes = c("TIRTLseqDataSet", "TIRTLseqData"))

  if(chain == "both") {
    df_a = forward_args(get_n_paired_step, chain = "alpha")
    df_b = forward_args(get_n_paired_step, chain = "beta")
    df_gg = bind_rows(df_a, df_b)
  } else {
    df_gg = forward_args(get_n_paired_step)
  }
  forward_args(make_step_plot_all, df_gg = df_gg)
}

#' Get a dataframe used for making the step plot (for all samples)
#' @noRd
get_n_paired_step = function(data, chain, cumulative = TRUE, freq = FALSE, step = "vh") {
  if("TIRTLseqDataSet" %in% class(data)) {
    df_all = lapply(names(data$data), function(x) {
      forward_args(get_n_paired_step_individual, data = data$data[[x]], chain = chain) |>
        mutate(sample = x, chain = chain)
    }) |> bind_rows()
  } else { ## TIRTLseqData
    df_all = forward_args(get_n_paired_step_individual, data = data, chain = chain) |>
      mutate(chain = chain)
  }
  return(df_all)
}

#' Make a ggplot with a step plot of number of partners for all samples
#' @noRd
make_step_plot_all = function(df_gg, chain, cumulative = TRUE, freq = FALSE, step = "vh", max_x = Inf) {
  if(!isTRUE(cumulative) && !isTRUE(freq)) ycol = "n_pairs_frac"; ylabel = "Number of pairs (fraction)"
  if(isTRUE(cumulative) && !isTRUE(freq)) ycol = "frac_cumulative"; ylabel = "Cumulative number of pairs (fraction)"
  if(!isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs"; ylabel = "Number of pairs"
  if(isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs_cumulative"; ylabel = "Cumulative number of pairs"
  xlabel = glue("Number of partners ({chain} chain)")
  max_x = min(max_x, max(df_gg$n))
  gg = ggplot(df_gg, aes(x=n, y=.data[[ycol]])) +
    geom_step(direction = step) +
    #scale_x_continuous( breaks = function(lims) c(1:5, seq(10, max(10, lims[2]), by = 5)), limits = c(0, max_x) ) +
    scale_x_continuous( breaks = function(lims) c(1:5, seq(10, max(10, lims[2]), by = 5)) ) +
    coord_cartesian(xlim = c(0, max_x)) +
    xlab(xlabel) + ylab(ylabel) +
    theme_bw()
  if ("sample" %in% colnames(df_gg)) gg = gg + aes(color = sample)
  if (chain == "both") gg = gg + facet_wrap(~chain)

  return(gg)
}

#' Get number of partners data frame for one sample
#' @noRd
get_n_paired_step_individual = function(data, chain = c("beta", "alpha"), cumulative = TRUE, freq = FALSE, step = "vh") {
  chain = chain[1]
  df_orig = data$paired_alt
  df = annotate_n_partners(df_orig)

  tab = table(df[[glue::glue("n_partners_{chain}")]])
  df_gg = tibble( n = as.integer(names(tab)), n_pairs = as.vector(tab) ) |>
    mutate(n_pairs_frac = n_pairs/sum(n_pairs)) |>
    mutate(frac_cumulative = cumsum(n_pairs_frac)) |>
    mutate(n_pairs_cumulative = cumsum(n_pairs)) |>
    bind_rows(data.frame(n=0L, n_pairs = 0L, n_pairs_frac = 0, frac_cumulative = 0, n_pairs_cumulative = 0L)) |>
    arrange(n)

  return(df_gg)
}



# plot_n_paired_step_individual = function(data, chain = c("beta", "alpha"), cumulative = TRUE, freq = FALSE, step = "vh") {
#   chain = chain[1]
#   df_orig = data$paired_alt
#   df = annotate_n_partners(df_orig)
#
#   tab = table(df[[glue::glue("n_partners_{chain}")]])
#   df_gg = tibble( n = as.integer(names(tab)), n_pairs = as.vector(tab) ) |>
#     mutate(n_pairs_frac = n_pairs/sum(n_pairs)) |>
#     mutate(frac_cumulative = cumsum(n_pairs_frac)) |>
#     mutate(n_pairs_cumulative = cumsum(n_pairs)) |>
#     bind_rows(data.frame(n=0L, n_pairs = 0L, n_pairs_frac = 0, frac_cumulative = 0, n_pairs_cumulative = 0L)) |>
#     arrange(n)
#
#   if(!isTRUE(cumulative) && !isTRUE(freq)) ycol = "n_pairs_frac"
#   if(isTRUE(cumulative) && !isTRUE(freq)) ycol = "frac_cumulative"
#   if(!isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs"
#   if(isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs_cumulative"
#   gg = ggplot(df_gg, aes(x=n, y=.data[[ycol]])) +
#     geom_step(direction = step) +
#     scale_x_continuous( breaks = function(lims) c(1:5, seq(10, max(10, lims[2]), by = 5)) ) +
#     theme_bw()
#
#   return(gg)
# }



# make_step_plot = function(df_gg, cumulative = TRUE, freq = FALSE, step = "vh", max_x = Inf) {
#   if(!isTRUE(cumulative) && !isTRUE(freq)) ycol = "n_pairs_frac"
#   if(isTRUE(cumulative) && !isTRUE(freq)) ycol = "frac_cumulative"
#   if(!isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs"
#   if(isTRUE(cumulative) && isTRUE(freq)) ycol = "n_pairs_cumulative"
#   max_x = min(max_x, max(df_gg$n))
#   gg = ggplot(df_gg, aes(x=n, y=.data[[ycol]], color = sample)) +
#     geom_step(direction = step) +
#     scale_x_continuous( breaks = function(lims) c(1:5, seq(10, max(10, lims[2]), by = 5)) ) +
#     theme_bw()
#   if ("sample" %in% colnames(df_gg)) gg <- gg + aes(color = sample)
#
#   return(gg)
# }





