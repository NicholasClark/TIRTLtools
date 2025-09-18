#' Scatterplot of TCR clone read fraction of clones between two samples
#'
#' @description
#' \code{plot_sample_vs_sample()} returns a scatterplot of read frequencies of TCRs between two samples
#'
#' @details
#' The function labels each TCR as up-regulated, down-regulated, or stable, based on the
#' log2 fold-change cutoff supplied (default 1.5).
#'
#' @param data1 a list of three data frames (alpha, beta, and paired) for one sample
#' @param data2 a list of three data frames (alpha, beta, and paired) for one sample
#' @param log2_cutoff the log2 fold-change cutoff to call a TCR up- or down-regulated (default 1.5)
#'
#' @return
#' A scatterplot (ggplot object) with read frequencies (proportions), colored by whether each
#' TCR is up-regulated, down-regulated, or neither, given the log2 fold-change cutoff.
#'
#' @family plotting
#' @export
#' @examples
#' # example code
#'
#'

plot_sample_vs_sample = function(data1, data2,
                                 chain = c("beta", "alpha"),
                                 #type_column = "auto",
                                 #value_type = c("auto","readFraction", "readCount", "n_wells", "readCount_max", "readCount_median", "avg", "n"),
                                 log2fc_cutoff = 3,
                                 sem_cutoff = 2.5,
                                 smooth_sem = c("window", "none"),
                                 window_size = 30,
                                 end_window_size = 5,
                                 pseudo1 = 1e-6,
                                 pseudo2 = 1e-6,
                                 labelx="Frequency on timepoint 1",
                                 labely="Frequency on timepoint 2",
                                 return_data = FALSE
                                 ) {
  chain = chain[1]
  smooth_sem = smooth_sem[1]
  # if(chain == "paired") {
  #   df1 = .add_single_chain_data_simple(data1) %>%
  #     mutate(readFraction = beta_readFraction,
  #            sem = beta_sem)
  #   df2 = .add_single_chain_data_simple(data2) %>%
  #     mutate(readFraction = beta_readFraction,
  #            sem = beta_sem)
  #   msg = paste("Using beta chain readFraction for frequencies") %>% .add_newline()
  #   cat(msg)
  # } else {
  #   df1 = data1[[chain]]
  #   df2 = data2[[chain]]
  # }

  df1 = data1[[chain]]
  df2 = data2[[chain]]

  # gg_df = .merge_TIRTL_pseudo(df1, df2, pseudo1 = pseudo1, pseudo2 = pseudo2)
  # #gg_df = .compute_fold_change(df1, df2, log2_cutoff = log2_cutoff, sem_cutoff = sem_cutoff)
  # gg_df = .add_sign(gg_df, pseudo1 = pseudo1, pseudo2 = pseudo2,
  #                  sem_threshold = sem_cutoff, log2FC_threshold = log2fc_cutoff
  #                  )
  if(smooth_sem == "window") {
    df1 = .smooth_sem_window(df1)
    df2 = .smooth_sem_window(df2)
  }
  gg_df = .compute_log2fc(df1, df2, pseudo1 = pseudo1, pseudo2 = pseudo2,
                         sem_cutoff = sem_cutoff, log2fc_cutoff = log2fc_cutoff
                         )
  gg = .plot_timepoints(gg_df, pseudo1 = pseudo1, pseudo2 = pseudo2,
                        labelx = labelx, labely = labely)
  if(return_data) return(gg_df)
  return(gg)
}

.compute_log2fc = function(df1, df2, pseudo1 = 1e-6, pseudo2 = 1e-6,
                          sem_cutoff = 2.5, log2fc_cutoff = 3
                          ) {
  df = .merge_TIRTL_pseudo(df1, df2, pseudo1 = pseudo1, pseudo2 = pseudo2)
  df = .add_sign(df, pseudo1 = pseudo1, pseudo2 = pseudo2,
                   sem_threshold = sem_cutoff, log2FC_threshold = log2fc_cutoff
  )
  return(df)
}

.moving_avg_closest <- function(x, window_size = 30, end_window_size=5) {
  n <- length(x)
  half_win <- floor(window_size / 2)


  out = sapply(1:length(x), function(i) {
    if(i >= length(x)-10) {
      half_win_end = floor(end_window_size / 2) ## smaller window for largest clones
      start = max(1, i - half_win_end)
      end   = min(n, i + half_win_end)
    } else {
      start = max(1, i - half_win)
      end   = min(n, i + half_win)
    }
    mean(x[start:end])
  })
}

.smooth_sem_window = function(df, window_size = 30, end_window_size=5) {
  df$sem_orig = df$sem
  df$sem = .moving_avg_closest(df$sem, window_size = window_size, end_window_size = end_window_size)
  #df$sem_log10 = .moving_avg_closest(log10(df$sem), window_size = window_size, end_window_size = end_window_size)
  #df$sem = 10^(df$sem_log10)
  return(df)
}

## df is the pseudobulk data frame (beta or alpha) for a sample
.plot_SEM_vs_read_fraction = function(df, window_size = 30, end_window_size=5) {
  if(!"sem_smoothed" %in% colnames(df)) df = .smooth_sem_window(df, window_size = window_size, end_window_size = end_window_size)

  ggplot(df, aes(x=log10(readFraction), y=log10(sem_orig))) +
    #geom_density_2d() +
    geom_hex(bins = 100) +
    geom_point(aes(x=log10(readFraction), y=log10(sem)), shape=1, alpha =0.5, color = "green") +
    #geom_point(shape = 1, alpha = 0.25) +
    #scale_x_log10() + scale_y_log10() +
    scale_fill_continuous(type = "viridis", limits = c(0, 500)) +
    theme_bw()
}

.plot_timepoints = function(dt,
                           labelx="Frequency on timepoint 1",
                           labely="Frequency on timepoint 2",
                           pseudo1 = 1e-6,
                           pseudo2 = 1e-6
                           ) {
  TIRTL_pallette = .get_tirtl_pallette()
  #col="#C5CBD3"
  dt$sign = factor(dt$sign, levels = c("stable", "down", "up"))
  log_labs_x = .get_log_labels_neg(dt$avg.x, pseudo1)
  log_labs_y = .get_log_labels_neg(dt$avg.y, pseudo2)
  ggplot(dt[dt$sign=="stable",], aes ((avg.x+pseudo1), (avg.y+pseudo2), color = sign))+
    geom_point(alpha=0.6, size=2)+
    geom_point(data=dt[dt$sign!="stable",], aes ((avg.x+pseudo1), (avg.y+pseudo2),color=sign), alpha=0.6, size=2)+
    theme_classic()+
    xlab(labelx)+
    ylab(labely)+
    geom_abline(linetype="dashed", col="black")+
    #  scale_color_manual(values=clrs[c(8,1,3)])+
    #theme(legend.position = "none")+
    scale_color_manual(values = c(stable = "grey70", down = TIRTL_pallette[10], up = TIRTL_pallette[7]))+
    scale_x_log10(breaks = log_labs_x$brks, labels = log_labs_x$labels) +
    scale_y_log10(breaks = log_labs_y$brks, labels = log_labs_y$labels)
    # scale_y_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),labels=c(expression("0"), expression("10"^"-5"), expression("10"^"-4"), expression("10"^"-3"), expression("10"^"-2")))+
    # scale_x_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),labels=c(expression("0"), expression("10"^"-5"), expression("10"^"-4"), expression("10"^"-3"), expression("10"^"-2")))#+
}

.merge_TIRTL_pseudo = function(tp1,tp2,thres1=4,thres2=4,pseudo1=1e-6,pseudo2=1e-6,mreads_thres=0)
{
  tmpm<-.na_to0(merge(tp1[readCount_median>mreads_thres,avg:=readFraction,],tp2[readCount_median>mreads_thres,avg:=readFraction,],by="targetSequences",all=T))[n_wells.x>thres1|n_wells.y>thres2,]
  tmpm[(n_wells.x<3),]$sem.x=mean(tmpm[n_wells.x==3,]$sem.x)*2
  tmpm[(n_wells.y<3),]$sem.y=mean(tmpm[n_wells.y==3,]$sem.y)*2
  tmpm[,log2FC:=log2((avg.y+pseudo1)/(avg.x+pseudo2)),]
  return(tmpm)
}

.add_sign<-function(tirtl_m,sem_threshold=2.5,log2FC_threshold=3,pseudo1=1e-6,pseudo2=1e-6){
  tirtl_m[,sign:="stable",]
  tirtl_m[log2FC<(-log2FC_threshold)&((avg.y+pseudo2+sem_threshold*sem.y)<(avg.x+pseudo1-sem_threshold*sem.x)),sign:="down",]
  tirtl_m[log2FC>(log2FC_threshold)&((avg.x+pseudo1+sem_threshold*sem.x)<(avg.y+pseudo2-sem_threshold*sem.y)),sign:="up",]
  return(tirtl_m)
}



# .compute_fold_change = function(data1, data2,
#                                 log2_cutoff = 3,
#                                 sem_cutoff = 5
# ) {
#   #value_type = value_type[1]
#   is_paired1 = .is.paired(data1)
#   is_paired2 = .is.paired(data2)
#   is_paired = is_paired1
#   if(is_paired1 != is_paired2) stop("Error: 'data1' and 'data2' need to be the same type -- either both paired or both single-chain.")
#   if(is_paired) join_col = "alpha_beta"
#   if(!is_paired) join_col = "targetSequences"
#
#   if(is_paired) {
#     data1 = remove_duplicates(data1)
#     data2 = remove_duplicates(data2)
#   }
#
#   join_df = inner_join(data1, data2, by = join_col)
#
#   proportion_column = "readFraction"
#   px = paste(proportion_column, ".x", sep = "")
#   py = paste(proportion_column, ".y", sep = "")
#   join_df$log2FC = log2(join_df[[py]]/join_df[[px]])
#   join_df$sem = pmax(join_df$sem.x, join_df$sem.y) ## take max of two SEMs?
#   join_df = join_df %>% mutate(
#     log2FC_signif = abs(log2FC) >= log2_cutoff,
#     sem_signif = abs(readFraction.x - readFraction.y) > 5*sem_cutoff
#     ) %>%
#     mutate(
#       sign = case_when(
#         (log2FC_signif & sem_signif) & (log2FC > 0) ~ "up",
#         (log2FC_signif & sem_signif) & (log2FC < 0) ~ "down",
#         .default = "stable"
#       )
#     )
#   join_df[["avg.x"]] = join_df[[px]]
#   join_df[["avg.y"]] = join_df[[py]]
#   return(join_df)
# }
#
#
# .compute_fold_change_old = function(data1, data2,
#                                type_column = "auto",
#                                value_type = c("auto","readFraction", "readCount", "n_wells", "readCount_max", "readCount_median", "avg", "n"),
#                                log2_cutoff = 1.5
# ) {
#   value_type = value_type[1]
#   is_paired1 = .is.paired(data1)
#   is_paired2 = .is.paired(data2)
#   is_paired = is_paired1
#   if(is_paired1 != is_paired2) stop("Error: 'data1' and 'data2' need to be the same type -- either both paired or both single-chain.")
#
#   proportion_column = .get_proportion_column(value_type, is_paired)
#   type_column = .get_type_column(type_column, is_paired)
#
#   cols = strsplit(type_column, "\\+")[[1]]
#   sym_type_col = syms(cols)
#
#   prop_df1 = data1 %>% group_by(!!!sym_type_col) %>%
#     summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE),
#               sem =
#               )
#   prop_df1$prop = prop_df1[[proportion_column]]/sum(prop_df1[[proportion_column]])
#   prop_df2 = data2 %>% group_by(!!!sym_type_col) %>%
#     summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
#   prop_df2$prop = prop_df2[[proportion_column]]/sum(prop_df2[[proportion_column]])
#
#   join_df = inner_join(prop_df1, prop_df2, by = cols)
#
#   px = paste(proportion_column, ".x", sep = "")
#   py = paste(proportion_column, ".y", sep = "")
#   join_df$log2FC = log2(join_df[[px]]/join_df[[py]])
#
#   join_df$sign = case_when(
#     join_df$log2FC >= log2_cutoff ~ "up",
#     join_df$log2FC <= -log2_cutoff ~ "down",
#     .default = "stable"
#   )
#   join_df[["avg.x"]] = join_df[[px]]
#   join_df[["avg.y"]] = join_df[[py]]
#   return(join_df)
# }
