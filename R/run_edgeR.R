
run_edgeR_list = function(
    ll,
    sample1_name,
    sample2_name,
    pseudobulk = TRUE,
    log2FC_cutoff = 3,
    fdr_cutoff = 0.05,
    reads_min = 1,
    n_wells_cutoff = 10,
    readFraction_cutoff = 1e-5,
    n_groups = 3,
    robust = TRUE,
    show_numbers = TRUE
) {
  mat1 = ll$matrix_list[[1]]
  mat2 = ll$matrix_list[[2]]
  colnames(mat1) = ll$col_data_list[[1]]$well
  colnames(mat2) = ll$col_data_list[[2]]$well
  row_data = ll$row_data
  col_data1 = ll$col_data_list[[1]]
  col_data2 = ll$col_data_list[[2]]
  res = run_edgeR(
    mat1 = mat1,
    mat2 = mat2,
    pseudobulk = pseudobulk,
    n_groups = n_groups,
    row_data = row_data,
    col_data1 = col_data1,
    col_data2 = col_data2,
    sample1_name = sample1_name,
    sample2_name = sample2_name,
    readFraction_cutoff = readFraction_cutoff,
    log2FC_cutoff = log2FC_cutoff,
    fdr_cutoff = fdr_cutoff,
    reads_min = reads_min,
    n_wells_cutoff = n_wells_cutoff,
    robust = robust,
    show_numbers = show_numbers
  )
  return(res)
}

run_edgeR = function(
    mat1,
    mat2,
    pseudobulk = TRUE,
    n_groups = 3,
    row_data,
    col_data1,
    col_data2,
    sample1_name,
    sample2_name,
    readFraction_cutoff = 1e-4,
    log2FC_cutoff = 3,
    fdr_cutoff = 0.05,
    reads_min = 10,
    n_wells_cutoff = 10,
    robust = TRUE,
    show_numbers = TRUE
    ) {

  mat1_full = mat1
  mat2_full = mat2
  ### make pseudo-bulk matrices
  if(pseudobulk) {
    ll1 = make_pseudo_replicates(mat1_full, n_groups)
    ll2 = make_pseudo_replicates(mat2_full, n_groups)
    mat1 = ll1$mat
    mat2 = ll2$mat
  }

  n1 = ncol(mat1) ## number of pseudo-replicates for group 1
  n2 = ncol(mat2) ## number of pseudo-replicates for group 2

  n_zeros_mat1 = Matrix::rowSums(mat1_full == 0)
  n_zeros_mat2 = Matrix::rowSums(mat2_full == 0)
  n_nonzero_mat1 = ncol(mat1_full) - n_zeros_mat1
  n_nonzero_mat2 = ncol(mat2_full) - n_zeros_mat2

  mat = Matrix::cbind2(mat1, mat2)
  grps = c(rep(sample1_name, n1), rep(sample2_name, n2))
  lib_sizes = Matrix::colSums(mat) ### keep track of total read counts for each sample
  read_frac1 = (Matrix::rowSums(mat1)+1)/sum(mat1) ### total read counts for each TCR (group 1)
  read_frac2 = (Matrix::rowSums(mat2)+1)/sum(mat2) ### total read counts for each TCR (group 2)
  log2FC = log2(read_frac2/read_frac1)
  abs_log2FC = abs(log2FC)

  #n_zeros_mat = Matrix::rowSums(mat == 0)
  #n_nonzero_mat = ncol(mat) - n_zeros_mat

  #junk = Matrix::rowSums(mat[which(n_nonzero_mat == 1),])

  # test_df = tibble(n_zeros = n_zeros_mat1, readFraction = read_frac1, n_nonzero = ncol(mat1)-n_zeros_mat1 )
  # ggplot(test_df) + geom_point(aes(x=n_nonzero, y=readFraction), alpha = 0.4) +
  #   scale_y_log10()
  # ggplot(tbl1 %>% filter(readCount1 > 0)) + geom_point(aes(x=n_nonzero1, y=readCount1), alpha = 0.4) +
  #   scale_y_log10()

  #keep_rows = which(total_counts >= sum(lib_sizes)*readFraction_cutoff) ### drop TCRs with total readFraction less than cutoff
  #keep_rows = which( ( (read_frac1 >= readFraction_cutoff) | (read_frac2 >= readFraction_cutoff) ) & (abs_log2FC >= log2FC_cutoff) )
  #keep_rows = which( (read_frac1 >= readFraction_cutoff) | (read_frac2 >= readFraction_cutoff) )




  row_data = row_data %>%
    mutate(readFraction1 = read_frac1, readFraction2 = read_frac2,
           log2FC_manual = log2FC, n_nonzero1 = n_nonzero_mat1,
           n_nonzero2 = n_nonzero_mat2, readCount1 = Matrix::rowSums(mat1),
           readCount2= Matrix::rowSums(mat2))


  #n_wells_cutoff1 = floor( (pct_wells_cutoff/100)*( ncol(mat1_full) ) )
  #n_wells_cutoff2 = floor( (pct_wells_cutoff/100)*( ncol(mat2_full) ) )
  n_wells_cutoff1 = n_wells_cutoff
  n_wells_cutoff2 = n_wells_cutoff

  #n_wells_cutoff = 3
  #keep_rows = (n_nonzero_mat1 >= n_wells_cutoff) | (n_nonzero_mat2 >= n_wells_cutoff)
  #reads_min = 5
  #keep_n_wells = ( Matrix::rowSums(mat1_full >= reads_min) >= n_wells_cutoff1 ) | ( Matrix::rowSums(mat2_full >= reads_min) >= n_wells_cutoff2 )
  keep_n_wells = ( Matrix::rowSums(mat1_full >= reads_min) >= n_wells_cutoff1 ) | ( Matrix::rowSums(mat2_full >= reads_min) >= n_wells_cutoff2 )
  keep_read_frac = (row_data$readFraction1 >= readFraction_cutoff) | (row_data$readFraction2 >= readFraction_cutoff)

  keep_rows = keep_n_wells & keep_read_frac

  mat_sub = mat[keep_rows,]
  row_data_sub = row_data[keep_rows,]
  rownames(mat_sub) = row_data_sub$row_id

  colnames(mat) = paste(grps, colnames(mat))
  colnames(mat_sub) =  paste(grps, colnames(mat_sub))

  # dge = edgeR::DGEList(counts = mat_sub, group = grps, lib.size = lib_sizes) ### need to add lib.size here for each sample
  # dge = edgeR::calcNormFactors(dge)


  #keep = filterByExpr(dge, group = dn_meta$timepoint)
  #keep = filterByExpr(dge)
  #dge = dge[keep, , keep.lib.sizes=FALSE]
  #dge = dge[keep, , keep.lib.sizes=TRUE]
  #keep <- filterByExpr(dge, group = sample_group, min.count = 2, min.total.count = 5)
  #dge <- dge[keep, , keep.lib.sizes = FALSE]


  dge = edgeR::DGEList(counts = mat_sub, group = grps, lib.size = lib_sizes)

  #### testing zinbwave
  # zinb = zinbFit(mat, K=2, epsilon=1000)
  # zinbwv = zinbwave(zinb, fitted_model = zinb, K = 2, epsilon=1000,
  #                  observationalWeights = TRUE)
  # weights <- assay(zinbwv, "weights")
  # dge <- DGEList(assay(fluidigm_zinb))
  # dge <- calcNormFactors(dge)
  #
  # design <- model.matrix(~Biological_Condition, data = colData(fluidigm))
  # dge$weights <- weights
  # dge <- estimateDisp(dge, design)
  # fit <- glmFit(dge, design)
  #
  # lrt <- glmWeightedF(fit, coef = 3)
  # topTags(lrt)
  ######## end testing zinbwave

  #keep1 = Matrix::rowSums(dge$counts > 0) >= 3
  #keep <- edgeR::filterByExpr(dge, group = grps, min.count = 1, min.total.count = 0)
  #dge <- dge[keep, , keep.lib.sizes = TRUE]
  dge = edgeR::calcNormFactors(dge)

  limma::plotMD(edgeR::cpm(dge, log=TRUE), column=1)
  abline(h=0, col="red", lty=2, lwd=2)
  md_plot1 = recordPlot()

  ### w/ design
  facs = factor(grps, levels = unique(grps))
  design = model.matrix(~0 + facs)
  colnames(design) = levels(facs)
  dge = edgeR::estimateDisp(dge, design, robust = robust) ## use robust = TRUE ??
  #fit = edgeR::glmFit(dge, design)
  print("BCV common dispersion")
  print(dge$common.dispersion)
  plotBCV(dge)
  bcv_plot = recordPlot()
  fit = edgeR::glmQLFit(dge, design, robust = robust) ## use robust = TRUE ??
  plotQLDisp(fit)
  qldisp = recordPlot()

  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = group2 - group1, levels = unique(grps))
  contrast_str <- paste0(sample2_name, " - ", sample1_name)
  contrast_label <- paste0(sample2_name, "_minus_", sample1_name)

  # Evaluate makeContrasts with a dynamically created formula
  call_str = paste0("limma::makeContrasts(",contrast_label, " = ", contrast_str, ", levels = design)")
  contrast1 = eval(parse(text = call_str))
  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = SUBJ023_CD8_d29 - SUBJ023_CD8_d2, levels = unique(grps))
  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = group2 - group1, levels = unique(grps))
  #lrt1 = edgeR::glmLRT(fit, contrast = contrast1)
  lrt1 = edgeR::glmQLFTest(fit, contrast = contrast1)

  de_top1 = edgeR::topTags(lrt1, n = Inf)
  de_df1 = de_top1$table %>%
    tibble::rownames_to_column("row_id") %>%
    mutate(row_id = as.integer(row_id)) %>%
    as_tibble() %>%
    left_join(row_data_sub, by = "row_id")

  call_df1 = limma::decideTests(lrt1, adjust.method = "BH", p.value = fdr_cutoff, lfc = log2FC_cutoff) %>%
    as.data.frame() %>% tibble::rownames_to_column("row_id") %>%
    set_colnames(c("row_id", "call")) %>%
    mutate(row_id = as.integer(row_id)) %>%
    as_tibble() %>%
    mutate(call_char = case_when(call == 0 ~ "stable", call == 1 ~ "up", call == -1 ~ "down"))

  de_df1 = de_df1 %>% left_join(call_df1, by = "row_id")
  up_df1 = de_df1 %>% filter(call == 1) %>% arrange(FDR)
  down_df1 = de_df1 %>% filter(call == -1) %>% arrange(FDR)
  non_sig = de_df1 %>% filter(call == 0) %>% arrange(FDR)

  limma::plotMD(lrt1, status = call_df1$call_char)
  abline(h=0, col="red", lty=2, lwd=2)
  md_plot2 = recordPlot()

  ### make traditional pseudo-bulk plot
  tbl1 = row_data %>% #mutate(readFraction1 = read_frac1, readFraction2 = read_frac2) %>%
    left_join(call_df1, by = "row_id") %>%
    mutate(call_char = ifelse(is.na(call_char), "not_tested", call_char))
  tbl1$call_char = factor(tbl1$call_char, levels = c("stable", "not_tested", "down", "up"))

  pseudo1 = 1e-6
  pseudo2 = 1e-6
  labelx = sample1_name
  labely = sample2_name
  TIRTL_pallette = .get_tirtl_pallette()

  n_up = nrow(up_df1)
  n_down = nrow(down_df1)
  n_tested = nrow(up_df1) + nrow(down_df1) + nrow(non_sig)

  log_labs_x = .get_log_labels_neg(tbl1$readFraction1, pseudo1)
  log_labs_y = .get_log_labels_neg(tbl1$readFraction2, pseudo2)
  gg = ggplot(tbl1[!tbl1$call_char %in% c("up", "down"),],
              aes ((readFraction1+pseudo1), (readFraction2+pseudo2), color = call_char))+
    geom_point(alpha=0.4, size=2)+
    geom_point(data=tbl1[tbl1$call_char %in% c("up", "down"),], aes ((readFraction1+pseudo1), (readFraction2+pseudo2),color=call_char), alpha=0.6, size=2)+
    theme_classic()+
    xlab(labelx)+
    ylab(labely)+
    geom_abline(linetype="dashed", col="black")+
    #  scale_color_manual(values=clrs[c(8,1,3)])+
    #theme(legend.position = "none")+
    scale_color_manual(values = c(stable = "grey70", not_tested = "grey60", down = TIRTL_pallette[10], up = TIRTL_pallette[7]))+
    scale_x_log10(breaks = log_labs_x$brks, labels = log_labs_x$labels) +
    scale_y_log10(breaks = log_labs_y$brks, labels = log_labs_y$labels)
  if(show_numbers) gg = gg + ggtitle(paste("Up:", n_up, "Down:", n_down, "Tested:", n_tested))
  print(gg)
  # gg = ggplot(tbl1) + geom_point(aes(x=readFraction1, y=readFraction2, color = call_char), alpha = 0.4) +
  #   geom_abline(slope = 1, intercept = 0) +
  #   scale_x_log10() +
  #   scale_y_log10() +
  #   theme_classic() +
  #   scale_color_manual(values = c(
  #     down = "steelblue",
  #     up = "firebrick",
  #     stable = "green",
  #     not_tested = "gray70"
  #   )) +
  #   xlab(paste(sample1_name, "read fraction")) +
  #   ylab(paste(sample2_name, "read fraction"))
  # print(gg)

  print(table(tbl1$call_char))

  percent_up_nonzero = (nrow(up_df1) - sum(up_df1$readFraction1 == min(de_df1$readFraction1, na.rm = TRUE) ))/nrow(up_df1)
  percent_down_nonzero = (nrow(down_df1) - sum(down_df1$readFraction2 == min(de_df1$readFraction2, na.rm = TRUE)))/nrow(down_df1)
  print("percent_up_nonzero")
  percent_up_nonzero %>% scales::label_percent(0.01)(.) %>% print()
  print("percent_up_nonzero")
  percent_down_nonzero %>% scales::label_percent(0.01)(.) %>% print()

  out_plots = list(md1 = md_plot1, md2 = md_plot2, bcv = bcv_plot, qldisp = qldisp, gg = gg)
  out_dfs = list(up = up_df1, down = down_df1, all = tbl1)
  out_vals = list(bcv_common_disp = dge$common.dispersion)

  res = list(plots = out_plots, dfs = out_dfs, values = out_vals)
  return(res)
}

## Split N wells into M groups
split_wells <- function(N, M) {
  base <- N %/% M        # minimum wells per group
  remainder <- N %% M    # number of groups that get an extra well
  rep(base + 1, remainder) |> c(rep(base, M - remainder))
}

## Make n pseudo-replicates from
make_pseudo_replicates = function(mat, n_groups) {
  wells_per_grp = split_wells(ncol(mat), n_groups)
  rand = sample(ncol(mat))
  grp_cols = split(rand, rep(1:n_groups, wells_per_grp))
  summed_groups <- sapply(grp_cols, function(cols) {
    Matrix::rowSums(mat[,cols, drop = F])
  })
  mat_new = Matrix::Matrix(summed_groups, sparse = TRUE)
  well_groups = lapply(grp_cols, function(i) colnames(mat)[i])
  return(list(mat = mat_new, groups = well_groups))
}
