

plot_timepoints = function(dt,
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
