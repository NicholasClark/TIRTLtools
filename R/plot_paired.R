### df is output from get_pair_stats()
plot_paired = function(df, chain = c("paired", "alpha", "beta")) {
  chain1 = chain[1]
  df_sub = df %>% filter(chain == chain1)
  ggplot(df_sub) + geom_col(aes(x=sample_id, y=Freq, fill = category))
}
