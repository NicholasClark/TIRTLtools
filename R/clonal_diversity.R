### clonal diversity

clonal_diversity = function(data, method = c("simpson","shannon")) {

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
