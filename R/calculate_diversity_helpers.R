### helper function - calculates diversity metrics for a single data frame
.diversity_single = function(data, chain, type_column = "auto", proportion_column="auto",
                             q=0:6, percent = seq(10,90,10), n = 10, tol = 1e-10,
                             methods = get_all_div_metrics()
) {
  data = data[[chain]]
  prop_df = .calculate_proportions(data = data, type_column = type_column, proportion_column = proportion_column)
  res = .calc_all_diversity(prop_df$prop, q=q, percent = percent, tol = tol, methods = methods, n = n)
  out = list(diversity = res, prop_df = prop_df)
  return(out)
}

### helper function to summarize clones by some column(s) and calculate proportions for each sample
# .calculate_proportions_list = function(data_list, type_column = "auto", proportion_column="auto", return_list = FALSE) {
#   is_paired = .is.paired(data_list)
#   if(type_column == "auto") {
#     if(is_paired) {
#       type_column = "alpha_beta"
#     } else {
#       type_column = "targetSequences"
#     }
#     msg = paste("\n", "Using ", type_column ," for 'type_column'", sep = "")
#     cat(msg)
#   }
#
#   cols = strsplit(type_column, "\\+")[[1]]
#   sym_type_col = syms(cols)
#
#   if(proportion_column == "auto") {
#     if(is_paired) {
#       proportion_column = "wij"
#     } else {
#       proportion_column = "readFraction"
#     }
#     msg = paste("\n", "Using ", proportion_column ," for 'proportion_column'", sep = "")
#     cat(msg)
#   }
#   print(length(data_list))
#   list_all = lapply(1:length(data_list), function(i) {
#     x = data_list[[i]]
#     tmp = .calculate_proportions(data = x, type_column = type_column, proportion_column = proportion_column) %>%
#       mutate(sample_num = i)
#     })
#   if(return_list) {
#     out = list_all
#   } else {
#     out = list_all %>% bind_rows()
#   }
#   return(out)
# }

### helper function to summarize clones by some column(s) for a single data frame
.calculate_proportions = function(data, type_column = "auto", proportion_column="auto") {
  #.calculate_proportions = function(data, type_column, proportion_column) {
  is_paired = .is.paired(data)
  if(.is.list.only(data)) stop("Error: 'data' needs to be a single data frame.")
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

  if(proportion_column == "auto") {
    if(is_paired) {
      proportion_column = "wij"
    } else {
      proportion_column = "readFraction"
    }
    msg = paste("\n", "Using ", proportion_column ," for 'proportion_column'", sep = "")
    cat(msg)
  }

  if(!proportion_column %in% colnames(data)) {
    warning("'proportion_column' not found in the data, using proportion of occurrences to measure diversity")
    #prop_df = data %>% group_by(!!sym(type_column)) %>% summarize(n = n())
    prop_df = data %>% group_by(!!!sym_type_col) %>% summarize(n = n())
    prop_df$prop = prop_df$n/sum(prop_df$n)
  } else {
    # prop_df = data %>% group_by(!!sym(type_column)) %>%
    #   summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
    prop_df = data %>% group_by(!!!sym_type_col) %>%
      summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
    prop_df$prop = prop_df[[proportion_column]]/sum(prop_df[[proportion_column]])
  }
  return(prop_df)
}


### helper function - calculates diversity metrics given a vector of proportions summing to one.
.calc_all_diversity = function(proportions, q=0:6, percent = seq(10,90,10), tol = 1e-14,
                               n=10,
                               methods = get_all_div_metrics() ) {
  data = proportions[!is.na(proportions)]
  data = data[data != 0]
  sum_data = sum(data)
  check1 = abs(sum_data - 1) < tol
  check2 = min(data) > 0 && max(data) <= 1
  if( (!check1) || (!check2) ) {
    s1 = "Data should be non-zero probabilities that sum to one."
    s2 = ""
    if(!check1) s2 = paste("Sum of data is:", sum_data)
    if(!check2) s2 = paste("Range of data is:", min(data), "to", max(data))
    ss = paste(c(s1, s2), collapse = " ")
    stop(ss)
  }
  ### return diversity indices for the requested methods
  res = list(
    simpson = .simpson(data),
    gini = .gini(data),
    gini.simpson = .gini.simpson(data),
    inv.simpson = .inv.simpson(data),
    shannon = .shannon(data),
    berger.parker = .berger.parker(data),
    richness = .richness(data),
    d50 = list(.d50(data)),
    dXX = list( .dxx_multi(p=data, xx=percent) ),
    renyi = list( .renyi_multi(p=data, q=q) ),
    hill = list( .hill_multi(p=data, q=q) ),
    top10fraction = .top10fraction(data),
    top100fraction = .top100fraction(data),
    topNfraction = list( .topNfraction(data, n) )
  )
  res = res[methods]
  return(res)
}
