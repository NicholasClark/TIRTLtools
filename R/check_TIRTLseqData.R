## Check that sample names are in the same order in the metadata and the data
check_TIRTLseqData = function(data, error = TRUE) {
  names1 = data$meta$sample_id
  names2 = names(data$data)
  check = identical(names1, names2)
  if(!check) {
    msg = "Names do not match between <data>$meta$sample_id and names(<data>$data)"
    if(error) {
      stop(msg)
    } else {
      warning(msg)
    }
  }
  return(check)
}

## input is a list with slots "paired", "alpha", and "beta"
.check_has_paired_data = function(data, error = TRUE, warn = TRUE) {
  check = isTRUE(tryCatch(nrow(data$paired) >= 1, error = function(e) FALSE))
  if(error && !check) stop("Error: object does not have paired data")
  if(warn && !check) warning("Warning: object does not have paired data")
  return(check)
}

## input is a list with slots "paired", "alpha", and "beta"
.check_has_single_chain_data = function(data, error = TRUE, warn = TRUE) {
  check1 = isTRUE(tryCatch({ nrow(data$alpha) >= 1 }, error = function(e) FALSE))
  check2 = isTRUE(tryCatch({ nrow(data$beta) >= 1}, error = function(e) FALSE))
  check = check1 && check2
  chain1 = ifelse(!check1, "alpha", "")
  chain2 = ifelse(!check2, "beta", "")
  chains = paste(chain1, chain2)
  test = error && !check
  if(error && !check) stop(paste("Error: object does not have single-chain data:", chains))
  if(warn && !check) warning(paste("Warning: object does not have single-chain data:", chains))
  return(check)
}
