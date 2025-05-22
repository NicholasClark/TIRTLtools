### all clonal diversity metrics

get_all_div_metrics = function() {
  div_metrics = c("simpson","gini","gini.simpson","inv.simpson","shannon",
    "berger.parker", "richness", "d50", "dXX", "renyi", "hill")
  return(div_metrics)
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
