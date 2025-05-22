#' @title
#' Calculation of clonal diversity measures from TCR data
#'
#' @description
#' `diversity` returns a number of clonal diversity measures (e.g. Hill numbers,
#' Shannon and Simpson indices) calculated from supplied data.
#'
#' @details
#' This function calculates the following diversity indices:
#'
#' simpson - The Simpson diversity index (https://en.wikipedia.org/wiki/Diversity_index#Simpson_index).
#' This equals the probability that two entities taken at random from the dataset represent the same type.
#' It is the inverse of the Hill number of order 2.
#'
#' gini - The Gini index/coefficient... (https://en.wikipedia.org/wiki/Gini_coefficient).
#'
#' gini.simpson - The Gini-Simpson index (https://en.wikipedia.org/wiki/Diversity_index#Gini%E2%80%93Simpson_index). Equal to 1-Simpson_index. This equals the probability that two entities taken at random from the dataset represent different types.
#'
#' inv.simpson - The Inverse-Simpson index - i.e. the reciprocal of the Simpson index, which measures the effective
#' number of types when the weighted arithmetic mean is used to calculate diversity. It is equivalent to the
#' Hill number of order 2. (https://en.wikipedia.org/wiki/Diversity_index#Inverse_Simpson_index)
#'
#' shannon - The Shannon diversity index, also known as the Shannon-Wiener index or Shannon entropy.
#' It is equivalent to the Renyi entropy when 'q' = 1.
#' (https://en.wikipedia.org/wiki/Diversity_index#Shannon_index)
#'
#' berger.parker - The Berger-Parker index, i.e. the proportion of the most abundant type in the dataset. (https://en.wikipedia.org/wiki/Diversity_index#Berger%E2%80%93Parker_index)
#'
#' richness - The species richness, i.e. the total number of unique types observed in the data (https://en.wikipedia.org/wiki/Diversity_index#Richness).
#'
#' d50 - The minimum number of types (clones) needed to comprise 50 percent of the data.
#'
#' dXX - The minimum number of types (clones) needed to comprise XX percent of the data.
#'
#' renyi - The Renyi entropy, a generalization of Shannon diversity/entropy for values of 'q' other than 1.
#' The formula for Renyi entropy is undefined at q=1, but it is defined as its limit, which equals the Shannon entropy.
#' When q=0, it is simply the natural logarithm of the richness or total number of types.
#' When q=Inf, it is defined by its limit, which is equal to the negative of the natural logarithm
#' of the proportion of the most abundant type.
#' (https://en.wikipedia.org/wiki/Diversity_index#R%C3%A9nyi_entropy)
#'
#' hill - The Hill numbers of order q, also known as the true diversity or the effective number of types.
#' This is the number of equally abundant types needed for the average proportional abundance of types to equal
#' that observed in the dataset. The order 'q' defines the exponent used in the formula to define the
#' 'generalized mean' of the proportional abundances.
#'
#' q = 2 corresponds to the arithmetic mean (https://en.wikipedia.org/wiki/Arithmetic_mean).
#'
#' q = 1 corresponds to the geometric mean (https://en.wikipedia.org/wiki/Geometric_mean).
#'
#' and q = 0 corresponds to the harmonic mean (https://en.wikipedia.org/wiki/Harmonic_mean).
#'
#' When q = Inf, it is defined by its limit, which is the inverse of the proportion of the most abundant type.
#'
#' In general, as q increases, less weight is given to rarer species.
#' (https://en.wikipedia.org/wiki/Diversity_index#Effective_number_of_species_or_Hill_numbers)
#'
#' @return returns a list with the requested indices. "simpson", "gini", "gini.simpson", "inv.simpson",
#' "shannon", "berger.parker", and "richness" return a vector of length one.
#'
#' "d50" and "dXX" return a data frame with the minimum number of clones needed to make up 50 (or XX) percent
#' of the data, the percentage (50 or XX) supplied, and the actual percentage of the data made up by those clones.
#' "dXX" may return a data frame with many rows if it is supplied a vector of percentages.
#'
#' "renyi" and "hill" return a data frame with the two columns, where each row contains the value for the Hill
#' (or Renyi) number of order 'q', and the other contains the value of the corresponding 'q'.
#' They return a data frame with many rows if supplied with a vector of orders (e.g. 1:5).
#'
#' @param data either a single data frame containing paired TCRs from TIRTL-seq output or a list of
#' data frames for many experiments
#' @param meta a dataframe of metadata for all of the samples in data. The first column is
#' expected to have unique sample names.
#' @param type_column the column from the TIRTLseq paired or pseudobulk output that determines
#' clonotype uniqueness. The default is "auto", which will use single-chain nucleotide sequence for pseudobulk and paired chain nucleotide sequence for paired data.
#' for paired results.
#' @param proportion_column the column from the TIRTLseq paired or pseudobulk data that will be used to
#' calculate the total proportion of each clonotype. The default is "readFraction". If no suitable
#' column is given, the proportion of occurences of each clonotype in the data will be used to
#' calculate diversity.
#' @param q a vector of integers specifying which "orders" to calculate for Hill numbers and Renyi entropy.
#' @param percent a percentage (out of 100) or a vector of percentages to use when calculating dXX values,
#' i.e. the minimum number of clones needed to cover XX percent of the sample.
#' @param tol the tolerance used to check that the proportions of all clones sum to 1 (default 10^-10)
#' @param methods the diversity indices or methods to calculate.
#'
#' @examples
#' # example here

#diversity = function(data, meta= NULL, chain = c("paired", "alpha", "beta"),
diversity = function(data, chain = c("paired", "alpha", "beta"),
                     type_column = "auto", proportion_column="auto",
                     q=0:6, percent = seq(10,90,10), tol = 1e-10,
                     methods = .get_all_div_metrics()
                     ) {
  meta = data$meta
  data = data$data
  chain = chain[1]
  if(!chain %in% names(data[[1]])) {
    msg = paste("Error: 'data' object supplied does not include data for chain: ", chain, "\n", sep = "")
    stop(msg)
  } else {
    msg = paste("Calculating diversity for chain: ", chain, "\n", sep = "")
    cat(msg)
  }
  # is_data_frame = is.data.frame(data[[chain]])
  # is_list = is.list(data[[chain]]) && !is_data_frame
  # is_paired = is.paired(data[[chain]])
  is_paired = chain == "paired"

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

  call_args = as.list(match.call())[-1]
  call_args$meta = NULL
  call_args$type_column = type_column
  call_args$proportion_column = proportion_column
  call_args$chain = chain

  n_samples = length(data)
  #if(is_list) {
    res = lapply(1:length(data), function(i) {
      msg = paste("\n", "-- Calculating diversity indices for sample ", i, " of ", n_samples,".", sep = "")
      cat(msg)
      # x=data[[chain]][[i]]
      x=data[[i]]
      call_args$data = x
      do.call(.diversity_single, call_args)
    }) %>% setNames(names(data))
  #} else {
  #  res = do.call(.diversity_single, call_args)
  #}
  out = list(result = res, meta = meta, call_args = call_args)
  return(out)
}

### helper function - calculates diversity metrics for a single data frame
.diversity_single = function(data, chain, type_column = "auto", proportion_column="auto",
                             q=0:6, percent = seq(10,90,10), tol = 1e-10,
                             methods = .get_all_div_metrics()
) {
  data = data[[chain]]
  prop_df = calculate_proportions(data = data, type_column = type_column, proportion_column = proportion_column)
  res = .calc_all_diversity(prop_df$prop, q=q, percent = percent, tol = tol, methods = methods)
  out = list(diversity = res, prop_df = prop_df)
  return(out)
}

### helper function to summarize clones by some column(s) and calculate proportions for each sample
calculate_proportions_list = function(data_list, type_column = "auto", proportion_column="auto", return_list = FALSE) {
  is_paired = is.paired(data_list)
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
  print(length(data_list))
  list_all = lapply(1:length(data_list), function(i) {
    x = data_list[[i]]
    tmp = calculate_proportions(data = x, type_column = type_column, proportion_column = proportion_column) %>%
      mutate(sample_num = i)
    })
  if(return_list) {
    out = list_all
  } else {
    out = list_all %>% bind_rows()
  }
  return(out)
}

### helper function to summarize clones by some column(s) for a single data frame
calculate_proportions = function(data, type_column = "auto", proportion_column="auto") {
#calculate_proportions = function(data, type_column, proportion_column) {
  is_paired = is.paired(data)
  if(is.list.only(data)) stop("Error: 'data' needs to be a single data frame.")
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
                            methods = .get_all_div_metrics() ) {
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
    hill = list( .hill_multi(p=data, q=q) )
  )
  res = res[methods]
  return(res)
}
