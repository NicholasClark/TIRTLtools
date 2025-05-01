### clonal diversity

clonal_diversity = function(data, q=5, percent=50, tol = 1e-14, method = c("simpson","shannon", "gini-simpson", "inverse-simpson", "renyi", "hill", "gini", "d50", "dXX")) {
  data = data[!is.na(data)]
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
  ### return div indices for the requested methods
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

.renyi = function(p, q) {
  (1/(1-q))*log(sum(p^q))
}

.hill = function(p, q) {
  if(q==0) {
    return( length(p) )
  } else if(q==1) {
    return ( exp(-sum(p*log(p))) ) ### exp(shannon entropy)
  } else if(q==Inf) {
    return ( 1/max(p) )
  } else {
    return ( sum(p^q)^(1/(1-q)) )
  }
}

## return number of species required to equal or surpass 'xx' percent of sample
.dxx = function(p, xx) {
  run_sum = cumsum(sort(p, decreasing = TRUE))
  idx = min(which(run_sum >= (xx/100)))
  return(idx)
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
