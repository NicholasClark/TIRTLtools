#' @title
#' Scatterplot of TCR frequency between two samples
#'
#' @description
#' \code{plot_sample_vs_sample()} returns a scatterplot of read frequencies of TCRs between two samples
#'
#' @details
#' The function labels each TCR as up-regulated, down-regulated, or stable, based on the
#' log2 fold-change cutoff supplied (default 1.5).
#'
#' @param data1 a dataframe for the first sample
#' @param data2 a dataframe for the second sample
#' @param log2_cutoff the log2 fold-change cutoff to call a TCR up- or down-regulated (default 1.5)
#'
#' @return
#' A scatterplot (ggplot object) with read frequencies (proportions), colored by whether each
#' TCR is up-regulated, down-regulated, or neither, given the log2 fold-change cutoff.
#'
#' @export
#' @examples
#' # example code
#'
#'

plot_sample_vs_sample = function(data1, data2,
                                 type_column = "auto",
                                 value_type = c("auto","readFraction", "readCount", "n_wells", "readCount_max", "readCount_median", "avg", "n"),
                                 log2_cutoff = 1.5) {
  gg_df = .compute_fold_change(data1, data2, type_column=type_column, value_type = value_type, log2_cutoff = log2_cutoff)
  gg = .plot_timepoints(gg_df)
  return(gg)
}

.plot_timepoints = function(dt,
                           labelx="Frequency on timepoint 1",
                           labely="Frequency on timepoint 2"
                           ) {
  TIRTL_pallet = c("#1D5F8A","#E76F47","#A53828","#531E1E","#116D4D","#ECC30B","#DA6140","#BABAA0","#38BBB7","#509E6E","grey90")

  ggplot(dt[dt$sign=="stable",], aes ((avg.x+1e-6), (avg.y+1e-6),col="#C5CBD3"))+
    geom_point(alpha=0.6, size=2)+
    geom_point(data=dt[dt$sign!="stable",], aes ((avg.x+1e-6), (avg.y+1e-6),col=sign), alpha=0.6, size=2)+
    theme_classic()+
    xlab(labelx)+
    ylab(labely)+
    geom_abline(linetype="dashed", col="black")+
    #  scale_color_manual(values=clrs[c(8,1,3)])+
    theme(legend.position = "none")+
    scale_color_manual(values = c("grey70",TIRTL_pallet[10],TIRTL_pallet[7]))+
    scale_y_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),labels=c(expression("0"), expression("10"^"-5"), expression("10"^"-4"), expression("10"^"-3"), expression("10"^"-2")))+
    scale_x_log10(breaks=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),labels=c(expression("0"), expression("10"^"-5"), expression("10"^"-4"), expression("10"^"-3"), expression("10"^"-2")))#+
}

.compute_fold_change = function(data1, data2,
                               type_column = "auto",
                               value_type = c("auto","readFraction", "readCount", "n_wells", "readCount_max", "readCount_median", "avg", "n"),
                               log2_cutoff = 1.5
) {
  value_type = value_type[1]
  proportion_column = value_type
  is_paired1 = .is.paired(data1)
  is_paired2 = .is.paired(data2)
  is_paired = is_paired1
  if(is_paired1 != is_paired2) stop("Error: 'data1' and 'data2' need to be the same type -- either both paired or both single-chain.")

  if(proportion_column == "auto") {
    if(is_paired) {
      proportion_column = "wij"
    } else {
      proportion_column = "readFraction"
    }
    msg = paste("\n", "Using ", proportion_column ," for 'proportion_column'", sep = "")
    cat(msg)
  }

  if(type_column == "auto") {
    if(is_paired) {
      type_column = "alpha_beta"
    } else {
      type_column = "targetSequences"
    }
    msg = paste("\n", "Using ", type_column ," for 'type_column'", sep = "")
    cat(msg)
  }
  cols = strsplit(type_column, "\\+")[[1]]
  sym_type_col = syms(cols)

  prop_df1 = data1 %>% group_by(!!!sym_type_col) %>%
    summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
  prop_df1$prop = prop_df1[[proportion_column]]/sum(prop_df1[[proportion_column]])
  prop_df2 = data2 %>% group_by(!!!sym_type_col) %>%
    summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
  prop_df2$prop = prop_df2[[proportion_column]]/sum(prop_df2[[proportion_column]])

  join_df = inner_join(prop_df1, prop_df2, by = cols)

  px = paste(proportion_column, ".x", sep = "")
  py = paste(proportion_column, ".y", sep = "")
  join_df$log2FC = log2(join_df[[px]]/join_df[[py]])
  join_df$sign = case_when(
    join_df$log2FC >= log2_cutoff ~ "up",
    join_df$log2FC <= -log2_cutoff ~ "down",
    .default = "stable"
  )
  join_df[["avg.x"]] = join_df[[px]]
  join_df[["avg.y"]] = join_df[[py]]
  return(join_df)
}
