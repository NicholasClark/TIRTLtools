#' Find expanded or contracted clones
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' 
#' This function returns clones that expand or contract between two samples, based
#' on their single-chain pseudo-bulk frequencies (beta chain frequency is default).
#' It returns data frames containing the most expanded/contracted clones along with their
#' direction and log-fold-change. It returns both a single-chain data frame (α or β, default is β)
#' and a data frame containing all αβ TCR pairs corresponding to these chains.
#' 
#' @details
#' If you would like a scatterplot of this data, you may call \code{\link{plot_sample_vs_sample}()}
#' with the same arguments for \code{log2fc_cutoff} and \code{sem_cutoff}.
#'
#' To call expanded and contracted clonotypes from TIRTL-seq data, we calculated mean frequency and
#' standard error of the mean (SEM) for each TCRβ chain over all wells. We call clones significantly
#' expanded or contracted between time points if there is a log2 fold-change log2FC > 3 between average 
#' frequencies and the difference between average frequencies exceeds 5 SEM intervals. This matches the
#' analysis in Pogorelyy & Kirk et al. (2025).
#' 
#' Note: For each TCR, we actually calculate two SEMs, one for each timepoint/sample.
#' To calculate 5 SEM intervals, we multiply each SEM by 2.5 and sum them. The \code{sem_cutoff}
#' argument controls this value, which is why the default is 2.5.
#'
#' @param data1 a sample from a TIRTLseqDataSet object (e.g. \code{<object>$data$<sample_tp1>}) - used for before frequencies
#' @param data2 a sample from a TIRTLseqDataSet object (e.g. \code{<object>$data$<sample_tp2>}) - used for after frequencies
#' @param chain which chain to plot, alpha or beta (default is beta)
#' @param filter_pairs whether to keep only two alpha chains per unique beta chain. (default is FALSE, keep all pairs)
#' @param remove_nonfunctional whether to remove pairs with non-functional chains --
#' e.g. that have stop codons or frameshifts. (default is FALSE, keep all pairs)
#' @param log2fc_cutoff the log2 fold-change cutoff to call a TCR expanded or contracted (default 1.5)
#' @param sem_cutoff the standard-error of the mean (SEM) to use as a cutoff in calling
#' clones expanded or contracted (default is 2.5)
#' @param pseudo1 the pseudocount to add to read frequency of the first sample (default is `10^-6`).
#' @param pseudo2 the pseudocount to add to read frequency of the second sample (default is `10^-6`).
#' @param smooth_sem if "window", then SEM values for clones will be smoothed by comparing to
#' other clones within a window of similar frequencies. Otherwise, no smoothing. (default is "window")
#' @param window_size the number of similar clones to include within a window.
#' @param end_window_size the number of clones to include in a window at the ends (most and least frequent)
#'
#' @returns
#' A list with two slots (`expanded` and `contracted`). Each slot is a list with two dataframes (`paired` and `single_chain`)
#' \preformatted{
#' (list)
#' └───expanded (list)
#'     └───paired (dataframe of all expanded αβ TCR pairs)
#'     └───single_chain (dataframe of all expanded single-chains)
#' └───contracted (list)
#'     └───paired (dataframe of all contracted αβ TCR pairs)
#'     └───single_chain (dataframe of all contracted single-chains)
#' }
#' 
#' @references
#' Pogorelyy, M, Kirk, A, Adhikari, S et al. (2025).
#' "TIRTL-seq: deep, quantitative and affordable paired TCR repertoire sequencing."
#' *Nature Methods*,
#' 23, 56–64. \doi{10.1038/s41592-025-02907-9}
#'
#' @family longitudinal
#' @export
#' @examples
#'
#' load_example_data() ## loads minimal SJTRC dataset into "SJTRC_minimal" object
#' clones = get_expanded_clones(
#'   data1 = SJTRC_minimal$data$cd8_tp1_v2,
#'   data2 = SJTRC_minimal$data$cd8_tp2_v2,
#'   chain = "beta",
#'   remove_nonfunctional = FALSE,
#'   filter_pairs = FALSE)
#' ## expanded clones
#' clones$expanded$paired ## paired data frame
#' clones$expanded$single_chain ## single_chain
#' ## contracted clones
#' clones$contracted$paired ## paired data frame
#' clones$contracted$single_chain ## single_chain
#' 

## testing
# load_example_data()
# data1 = SJTRC_minimal$data$cd8_tp1_v2
# data2 = SJTRC_minimal$data$cd8_tp2_v2
# chain = "beta"
# filter_pairs = FALSE
# remove_nonfunctional = FALSE
# log2fc_cutoff = 3
# sem_cutoff = 2.5
# pseudo1 = 1e-6
# pseudo2 = 1e-6
# smooth_sem = c("window", "none")
# window_size = 30
# end_window_size = 5

# load_example_data()
# test = get_expanded_clones(SJTRC_minimal$data$cd8_tp1_v2, SJTRC_minimal$data$cd8_tp2_v2, remove_nonfunctional = F, filter_pairs = F)

get_expanded_clones = function( data1,
                                data2,
                                chain = c("beta", "alpha"),
                                filter_pairs = FALSE,
                                remove_nonfunctional = TRUE,
                                log2fc_cutoff = 3,
                                sem_cutoff = 2.5,
                                pseudo1 = 1e-6,
                                pseudo2 = 1e-6,
                                smooth_sem = c("window", "none"),
                                window_size = 30,
                                end_window_size = 5) {
  chain = chain[1]
  smooth_sem = smooth_sem[1]
  checkmate::assert_choice(chain, choices = c("alpha", "beta"))
  
  other_chain = ifelse(chain == "alpha", "beta", "alpha")

  df_all = .make_sample_vs_sample_df(data1 = data1, data2 = data2, chain = chain, smooth_sem = smooth_sem,
    pseudo1=pseudo1, pseudo2 = pseudo2, sem_cutoff = sem_cutoff, log2fc_cutoff = log2fc_cutoff,
    window_size = window_size, end_window_size = end_window_size
  )
  paired1 = data1$paired_alt %>% mutate(pairs_from = "x")
  paired2 = data2$paired_alt %>% mutate(pairs_from = "y")

  up_df = df_all %>% filter(sign %in% c("up")) %>%
    arrange(desc(abs(log2FC))) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())
  down_df = df_all %>% filter(sign %in% c("down")) %>%
    arrange(desc(abs(log2FC))) %>%
    select(targetSequences, aaSeqCDR3, v, j, sign, log2FC, avg.x, avg.y, sem.x, sem.y, everything())

  up_out = .merge_with_pairs(single_chain_df = up_df, direction = "up", paired_df1 = paired1, paired_df2 = paired2,
    chain = chain, filter_pairs = filter_pairs, remove_nonfunctional = remove_nonfunctional)
  down_out = .merge_with_pairs(single_chain_df = down_df, direction = "down", paired_df1 = paired1, paired_df2 = paired2,
    chain = chain, filter_pairs = filter_pairs, remove_nonfunctional = remove_nonfunctional)  

  out_list = list(expanded = up_out, contracted = down_out)
  return(out_list)
}

.merge_with_pairs = function(
  single_chain_df, direction, paired_df1, paired_df2,
  chain, filter_pairs, remove_nonfunctional) {
  
  if(chain == "alpha") {
    cdr3_nt_col = "alpha_nuc"
    cdr3_aa_col = "cdr3a"
    v_col = "va"
    j_col = "ja"
  } else {
    cdr3_nt_col = "beta_nuc"
    cdr3_aa_col = "cdr3b"
    v_col = "vb"
    j_col = "jb"
  }
  join_cols = c(cdr3_nt_col, cdr3_aa_col, v_col, j_col)

  df_select = single_chain_df %>% 
    select(
      sign, log2FC, 
      targetSequences, aaSeqCDR3, v, j,
      #avg.x, avg.y, 
      readFraction.x, readFraction.y,
      sem.x, sem.y, n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readCount.x, readCount.y,
      readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y
      ) %>%
    #dplyr::rename(beta_nuc = targetSequences, cdr3b = aaSeqCDR3, vb = v, jb = j)
    dplyr::rename(!!cdr3_nt_col := targetSequences, !!cdr3_aa_col := aaSeqCDR3, !!v_col := v, !!j_col := j)
  
  if(filter_pairs) {
    ## filter paired data frame from each sample down to the best partner for each beta
    paired_df1 = .clean_pairs_single(paired_df1)
    paired_df2 = .clean_pairs_single(paired_df2)
  }
  paired_df1_select = paired_df1 %>% 
    select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from)
  paired_df2_select = paired_df2 %>% 
    select(va, ja, cdr3a, cdr3b, vb, jb,
      alpha_nuc, beta_nuc, alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from)
  
  
  signif_paired_df1 = df_select %>% inner_join(paired_df1_select, by = join_cols)
  signif_paired_df2 = df_select %>% inner_join(paired_df2_select, by = join_cols)
  if(direction == "up") {
    ## if clones are expanding put pairs from second sample at the top (more likely to find pair there)
    signif_paired_all = bind_rows(signif_paired_df2, signif_paired_df1)
  } else {
    ## if clones are contracting put pairs from first sample at the top (more likely to find pair there)
    signif_paired_all = bind_rows(signif_paired_df1, signif_paired_df2)
  }

  signif_paired_all = signif_paired_all %>%
    mutate(comparison_chain = chain) %>%
    filter(!duplicated(alpha_beta)) %>%
    select(
      sign, log2FC,
      va, ja, cdr3a, cdr3b, vb, jb, #avg.x, avg.y,
      alpha_nuc, beta_nuc,
      comparison_chain,
      readFraction.x, readFraction.y,
      sem.x, sem.y, 
      n_wells.x, n_wells.y, max_wells.x, max_wells.y,
      rank.x, rank.y, readCount.x, readCount.y,
      readCount_median.x, readCount_median.y, readCount_max.x, readCount_max.y,
      alpha_beta, 
      is_functional, alpha_is_functional, beta_is_functional, pairs_from) %>%
    mutate(across(starts_with("readCount") | starts_with("max_wells") | 
      starts_with("n_wells") | starts_with("rank"), as.integer)) %>%
    arrange(desc(abs(log2FC))) %>%
    as_tibble()
  
  ## single-chain output
  signif_sc = df_select %>% 
    mutate(across(starts_with("readCount") | starts_with("max_wells") | 
      starts_with("n_wells") | starts_with("rank"), as.integer)) %>%
    mutate(is_paired = !!sym(cdr3_nt_col) %in% signif_paired_all[[cdr3_nt_col]]) %>%
    #filter(!(!!sym(cdr3_nt_col) %in% signif_paired_all[[cdr3_nt_col]])) %>% ## filter out paired single-chains
    select(sign, log2FC, is_paired, everything()) %>%
    arrange(desc(abs(log2FC))) %>%
    as_tibble()
  if(remove_nonfunctional) {
    signif_paired_all = signif_paired_all %>% filter(is_functional)
    signif_sc = signif_sc %>%
      mutate(has_stop_codon = grepl("\\*", !!sym(cdr3_aa_col)), has_frameshift = grepl("_", !!sym(cdr3_aa_col))) %>%
      mutate(is_functional = !(has_stop_codon | has_frameshift)) %>%
      filter(is_functional)
  }
  out = list(paired = signif_paired_all, single_chain = signif_sc)
  return(out)
}
