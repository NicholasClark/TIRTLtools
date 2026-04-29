#' Plot a plate of single-cell TIRTL-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This function plots which wells in a plate are missing an alpha and/or beta chain
#' for single-cell TIRTL-seq data.
#'
#' @param df a data frame, output from the \code{\link{process_scTIRTLseq}()} function
#' @param title an optional title for the plot
#'
#' @return
#' A ggplot object.
#' @family single-cell
#' @export
#'

plot_scTIRTLseq_plate = function(df, title = "") {
  TIRTL_pallet2<-c("#007EA7","#3CAF82","#FFD16E","#A53828","#F37748","grey90")
  substitution3 <- c("TRUETRUE" = "ab", "TRUEFALSE" = "a lost", "FALSETRUE" = "b lost")
  df2 = df[,.(barcode,chains=substitution3[paste0(!is.na(nSeqCDR3_beta),!is.na(nSeqCDR3_alpha))]),]
  gg = ggplate::plate_plot(df2,
                  plate_size = 384,
                  position = barcode,
                  value=chains,show_legend = FALSE,title=NULL) +
    scale_fill_manual(values=c("ab"=TIRTL_pallet2[2],"b lost"=TIRTL_pallet2[5],"a lost"=TIRTL_pallet2[3])) +
    theme(legend.position = "bottom") +
    ggtitle(title)

  return(gg)
}
