#' Returns all diversity metric options for \code{\link{calculate_diversity}()}
#'
#' @return
#' A vector of all available diversity metrics for the \code{\link{calculate_diversity}()} function.
#'
#' @family repertoire_analysis
#' @seealso \code{\link{diversity}()}, \code{\link{plot_diversity}()}
#'
#' @export
#' @examples
#' # example code
#' get_all_div_metrics()
#'
get_all_div_metrics = function() {
  div_metrics = c("simpson","gini","gini.simpson","inv.simpson","shannon",
                  "berger.parker", "richness", "d50", "dXX", "renyi", "hill", "top10fraction",
                  "top100fraction", "topNfraction")
  return(div_metrics)
}
