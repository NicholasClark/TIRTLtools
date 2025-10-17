
run_edgeR = function(mat1, mat2, row_data, col_data1, col_data2, sample1_name, sample2_name,
                     readFraction_cutoff = 1e-6, log2FC_cutoff = 2, fdr_cutoff = 0.05) {
  n1 = ncol(mat1)
  n2 = ncol(mat2)
  mat = Matrix::cbind2(mat1, mat2)
  grps = c(rep(sample1_name, n1), rep(sample2_name, n2))
  lib_sizes = Matrix::colSums(mat) ### keep track of total read counts for each sample
  read_frac1 = (Matrix::rowSums(mat1)+1)/sum(mat1) ### total read counts for each TCR (group 1)
  read_frac2 = (Matrix::rowSums(mat2)+1)/sum(mat2) ### total read counts for each TCR (group 2)
  log2FC = log2(read_frac2/read_frac1)
  abs_log2FC = abs(log2FC)
  #keep_rows = which(total_counts >= sum(lib_sizes)*readFraction_cutoff) ### drop TCRs with total readFraction less than cutoff
  keep_rows = which( ( (read_frac1 >= readFraction_cutoff) | (read_frac2 >= readFraction_cutoff) ) & (abs_log2FC >= log2FC_cutoff) )
  mat_sub = mat[keep_rows,]
  row_data_sub = row_data[keep_rows,]
  rownames(mat_sub) = row_data_sub$row_id

  dge = edgeR::DGEList(counts = mat_sub, group = grps, lib.size = lib_sizes) ### need to add lib.size here for each sample
  #keep = filterByExpr(dge, group = dn_meta$timepoint)
  #keep = filterByExpr(dge)
  #dge = dge[keep, , keep.lib.sizes=FALSE]
  #dge = dge[keep, , keep.lib.sizes=TRUE]
  dge = edgeR::calcNormFactors(dge)

  ### w/ design
  facs = factor(grps, levels = unique(grps))
  design = model.matrix(~0 + facs)
  colnames(design) = levels(facs)
  dge = edgeR::estimateDisp(dge, design)
  fit = edgeR::glmFit(dge, design)

  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = group2 - group1, levels = unique(grps))
  contrast_str <- paste0(sample2_name, " - ", sample1_name)
  contrast_label <- paste0(sample2_name, "_minus_", sample1_name)

  # Evaluate makeContrasts with a dynamically created formula
  call_str = paste0("limma::makeContrasts(",contrast_label, " = ", contrast_str, ", levels = design)")
  contrast1 = eval(parse(text = call_str))
  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = SUBJ023_CD8_d29 - SUBJ023_CD8_d2, levels = unique(grps))
  #contrast1 = limma::makeContrasts(grp2_minus_grp1 = group2 - group1, levels = unique(grps))
  lrt1 = edgeR::glmLRT(fit, contrast = contrast1)

  de_top1 = edgeR::topTags(lrt1, n = Inf)
  de_df1 = de_top1$table %>% tibble::rownames_to_column("row_id") %>%
    mutate(row_id = as.integer(row_id)) %>%
    as_tibble() %>%
    left_join(row_data_sub, by = "row_id")

  limma::plotMD(lrt1)
  md_plot = recordPlot()

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


  ### make traditional pseudo-bulk plot
  tbl1 = row_data %>% mutate(readFraction1 = read_frac1, readFraction2 = read_frac2) %>%
    left_join(call_df1, by = "row_id") %>%
    mutate(call_char = ifelse(is.na(call_char), "not_tested", call_char))

  gg = ggplot(tbl1) + geom_point(aes(x=readFraction1, y=readFraction2, color = call_char), alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    scale_color_manual(values = c(
      down = "steelblue",
      up = "firebrick",
      stable = "gray40",
      not_tested = "gray70"
    )) +
    xlab(paste(sample1_name, "read fraction")) +
    ylab(paste(sample2_name, "read fraction"))
  gg

  table(tbl1$call_char)

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
