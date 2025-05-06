### clonal diversity

#' possible additions:
#' 1) chao diversity (https://www.uvm.edu/~ngotelli/manuscriptpdfs/ChaoHill.pdf), rarefaction (https://en.wikipedia.org/wiki/Rarefaction_(ecology))
#' 2) plotting
#' 3) (done) allow to take in a list of data frames
#' 4) only calculate specified metrics to make it faster
#' 5) allow for multiple columns to be used as "type" -- i.e. to group by both alpha and beta chains (va, vb, cdr3a, cdr3b)



#' @title
#' Calculation of clonal diversity measures from TCR data
#'
#' @description
#' `clonal_diversity` returns a number of clonal diversity measures (e.g. Hill numbers,
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
#' gini.simpson - The Gini-Simpson index (https://en.wikipedia.org/wiki/Diversity_index#Gini%E2%80%93Simpson_index)
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
#' @param q a vector of integers specifying which "orders" to calculate for Hill numbers and Renyi entropy.
#' @param percent a percentage (out of 100) or a vector of percentages to use when calculating dXX values,
#' i.e. the minimum number of clones needed to cover XX percent of the sample.
#' @param tol the tolerance used to check that the proportions of all clones sum to 1
#' @param methods the diversity indices or methods to calculate.
#'
#' @examples
#'

diversity = function(data, type_column = "auto", proportion_column="readFraction",
                     q=0:6, percent = seq(10,90,10), tol = 1e-10,
                     methods = c("simpson","gini","gini.simpson","inv.simpson","shannon",
                                 "berger.parker", "richness", "d50", "dXX", "renyi", "hill")
                     ) {

  is_data_frame = is.data.frame(data)
  is_list = is.list(data) && !is_data_frame
  if(!(is_list || is_data_frame)) stop("'data' needs to be a data frame or a list of data frames")
  if(is_data_frame) is_paired = "wij" %in% colnames(data)
  if(is_list) is_paired = "wij" %in% colnames(data[[1]])



  call_args <- as.list(match.call())[-1]

  if(is_list) {
    res = lapply(data, function(x) {
      call_args$data = x
      #tmp = .diversity_single(x, type_column = type_column, )
      do.call(.diversity_single, call_args)
    })
  } else {
    res = do.call(.diversity_single, call_args)
  }
  return(res)
}

.diversity_single = function(data, type_column = "auto", proportion_column="readFraction",
                             q=0:6, percent = seq(10,90,10), tol = 1e-10,
                             methods = c("simpson","gini","gini.simpson","inv.simpson","shannon",
                                         "berger.parker", "richness", "d50", "dXX", "renyi", "hill")
) {
  if(type_column == "auto") {
    if(is_paired) {
      type_column = "cdr3b"
    } else {
      type_column = "aaSeqCDR3"
    }
  }
  if(!proportion_column %in% colnames(data)) {
    warning("'proportion_column' not found in the data, using proportion of occurrences to measure diversity")
    prop_df = data %>% group_by(!!sym(type_column)) %>% summarize(n = n())
    prop_df$prop = prop_df$n/sum(prop_df$n)
  } else {
    prop_df = data %>% group_by(!!sym(type_column)) %>%
      summarize(!!sym(proportion_column) := sum(!!sym(proportion_column), na.rm = TRUE))
    prop_df$prop = prop_df[[proportion_column]]/sum(prop_df[[proportion_column]])
  }
  res = .calc_all_diversity(prop_df$prop, q=q, percent = percent, tol = tol, methods = methods)
  return(res)
}

plot

.calc_all_diversity = function(proportions, q=0:6, percent = seq(10,90,10), tol = 1e-14,
                            methods = c("simpson","gini","gini.simpson","inv.simpson","shannon",
                                        "berger.parker", "richness", "d50", "dXX", "renyi", "hill")) {
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

.berger.parker = function(p) {
  max(p)
}

.richness = function(p) {
  length(p)
}

.simpson = function(p) {
  sum(p^2)
}

.gini.simpson = function(p) {
  1-sum(p^2) ### 1 - simpson
}

.inv.simpson = function(p) {
  1/sum(p^2) ### 1/simpson
}

.shannon = function(p) {
  -sum(p*log(p))
}

.renyi_multi = function(p,q) {
  res = lapply(q, function(x) .renyi(p,x))
  return(bind_rows(res))
}

.renyi = function(p, q) {
  if(length(q) > 1) stop("'q' needs to be a single number")
  if(q==1) return( list (value = .shannon(p), q = q) )
  if(q==Inf) return( list( value = -log(max(p)), q = q ) )
  return( list( value = (1/(1-q))*log(sum(p^q)), q = q) )
}

.hill_multi = function(p, q) {
  res = lapply(q, function(x) .hill(p,x))
  return(bind_rows(res))
}

.hill = function(p, q) {
  if(length(q) > 1) stop("'q' needs to be a single number")
  if(q==0) {
    return( list(value = length(p), q = q) )
  } else if(q==1) {
    return ( list(value = exp(-sum(p*log(p))), q = q) ) ### exp(shannon entropy)
  } else if(q==Inf) {
    return ( list(value = 1/max(p), q = q) )
  } else {
    return ( list(value = sum(p^q)^(1/(1-q)), q = q) )
  }
}

.dxx_multi = function(p, xx_vec) {
  res = lapply(xx_vec, function(x) .dxx(p,x))
  return(bind_rows(res))
}

## return number of species required to equal or surpass 'xx' percent of sample
.dxx = function(p, xx) {
  if(length(xx) > 1) stop("'xx' needs to be a single number")
  if(xx < 0 || xx > 100) stop("'xx' needs to be a percentage between 0 and 100")
  run_sum = cumsum(sort(p, decreasing = TRUE))
  idx = min(which(run_sum >= (xx/100)))
  return( tibble(n_clones = idx, percent_required = xx, percent_observed = 100*run_sum[idx]) )
}

.gini = function(p) {
  p_sort = sort(p, decreasing = FALSE)
  n = length(p)
  g = (1/n)*(n+1-2*(sum( (n:1)*p_sort )/sum(p_sort)))
  return(g)
}

## return number of species required to equal or surpass 50 percent of sample
.d50 = function(p) {
  .dxx(p, 50)
}

# rarefaction = function(p, step, quantile) {
#
# }
