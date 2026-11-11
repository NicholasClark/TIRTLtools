

get_96well_subset = function(subset = c("A1", "B1", "A2", "B2")) {
  subset = subset[1]
  checkmate::assert_choice(subset, c("A1", "B1", "A2", "B2"))
  row_start = substr(subset, 0, 1)
  row_start_int = if(row_start == "A") 1 else 2
  col_start_int = substr(subset, 2, 3) |> as.numeric()
  wells = unlist(sapply(LETTERS[seq(from=row_start_int,to=16,by=2)],function(x)paste(x,seq(from=col_start_int,to=24,by=2),sep=""),simplify = F)) |> as.vector()
  return(wells)
}
