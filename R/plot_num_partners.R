plot_num_partners = function(data,
                             label_col = NULL, group_col = NULL, fraction = TRUE,
                             include_non_functional = FALSE,
                             max_partners = 5, return_data = FALSE) {

  meta = data$meta
  data = data$data

  if(include_non_functional) {
    data = lapply(data, function(x) x[["paired"]]) %>% setNames(names(data))
  } else {
    data = lapply(data, function(x) x[["paired"]] %>% filter(is_functional)) %>% setNames(names(data))
  }


  is_paired = is.paired(data)
  is_list = is.list.only(data)
  if(!is_paired) stop("'data' must be paired chain output from TIRTL-seq")
  data = remove_dupes_paired(data)

  if(is_list) {
    gg_df = lapply(1:length(data), function(i) {
      df_tmp = get_num_partners_single(data[[i]], max_partners = max_partners)
      df_tmp$sample_num = i
      if(!is.null(group_col)) {
        df_tmp$Group = meta[[group_col]][i]
      } else {
        df_tmp$Group = meta[[1]][i]
      }
      if(!is.null(label_col)) df_tmp$label = meta[[label_col]][i]
      return(df_tmp)
      }) %>% bind_rows() %>% group_by(Group, n_partners, chain) %>%
      summarize(Frequency = sum(Frequency), Fraction = mean(Fraction)) ## note this is taking the mean proportion over samples rather than summing all frequencies and dividing by the total
  } else {
    gg_df = get_num_partners_single(data)
  }
  var = ifelse(fraction, sym("Fraction"), sym("Frequency"))
  char = paste(">", max_partners, sep = "")
  lvls = c(1:max_partners, char) %>% rev()
  gg_df$n_partners = factor(gg_df$n_partners, levels = lvls)
  gg = ggplot(gg_df) +
    geom_col(aes(x = chain, y = !!var, fill = n_partners)) +
    xlab("") +
    theme_classic()
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

get_num_partners_single = function(df, max_partners = 5) {
  alpha_tbl = table(df$alpha_nuc) %>% table() %>% as.data.frame.table() %>%
    magrittr::set_colnames(c("n_partners", "Frequency")) %>%
    mutate( Fraction = Frequency/sum(Frequency), chain = "alpha" )
  beta_tbl = table(df$beta_nuc) %>% table() %>% as.data.frame.table() %>%
    magrittr::set_colnames(c("n_partners", "Frequency")) %>%
    mutate( Fraction = Frequency/sum(Frequency), chain = "beta" )
  df_long = bind_rows(alpha_tbl, beta_tbl) %>% mutate(n_partners = as.integer(as.character(n_partners)))
  gt_char = paste(">", max_partners, sep = "")
  df_less_max = df_long %>% filter(n_partners <= max_partners) %>% mutate(n_partners = as.character(n_partners))
  df_greater = df_long %>% filter(n_partners > max_partners) %>%
    group_by(chain) %>%
    summarize(Frequency = sum(Frequency), Fraction = sum(Fraction)) %>%
    ungroup() %>%
    mutate(n_partners = gt_char)
  df_out = bind_rows(df_less_max, df_greater) #%>%
    #mutate(n_partners = factor(n_partners, levels = c(as.character(sort(as.integer(unique(df_less_max$n_partners)))), gt_char)) )
  return(df_out)
}
