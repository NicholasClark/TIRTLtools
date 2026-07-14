#' Load example TIRTL-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' 
#' This is a helper function to load example TIRTL-seq data. The function will also
#' download the data if needed from \link{https://github.com/NicholasClark/TIRTLtools/releases/tag/data-v1}
#' and save it in the package's R cache directory.
#' 
#' The function currently loads one of two available TIRTL-seq datasets, 
#' either "SJTRC_minimal" or "SJTRC_longitudinal". The data is loaded
#' to the global environment and assigned with the same name as the dataset.
#' For example, if \code{dataset} is "SJTRC_minimal", the data will be set to `SJTRC_minimal`
#' in the global environment. If the dataset is already loaded, it will not be re-loaded.
#'
#' @details
#' The \code{SJTRC_longitudinal} dataset contains 6 samples of TIRTL-seq data from the
#' St. Jude Tracking Study of Immune Responses Associated with COVID-19 (SJTRC).
#' The data is from one donor at three timepoints, isolated for either CD4+ or CD8+ T-cells.
#'
#' The \code{SJTRC_minimal} dataset is a subset of the above dataset, containing only two samples:
#' CD8+ isolated T-cells for timepoints 1 and 2.
#'
#' @param dataset either "SJTRC_minimal" or "SJTRC_longitudinal".
#' @param verbose whether to print messages to the user (default TRUE)
#'
#' @return
#' The function returns NULL invisibly, but loads a dataset and assigns it to either 
#' `SJTRC_minimal` or `SJTRC_longitudinal` in the global environment.
#'
#' @family data_loading
#'
#' @examples
#' load_example_data("SJTRC_minimal")
#' print(SJTRC_minimal)
#' summary(SJTRC_minimal)
#' load_example_data("SJTRC_longitudinal")
#' print(SJTRC_longitudinal)
#' summary(SJTRC_longitudinal)
#'

load_example_data = function(dataset = c("SJTRC_minimal", "SJTRC_longitudinal"), verbose = TRUE) {
  dataset = dataset[1]
  assert_choice(dataset, choices = c("SJTRC_minimal", "SJTRC_longitudinal"))
  assert_logical(verbose)
  if (!interactive()) verbose = FALSE

  if (exists(dataset, envir = .GlobalEnv, inherits = FALSE)) {
    if(verbose) message(glue("Example data already loaded into object: '{dataset}'"))
    return(invisible(NULL))
  }
  ts_data = load_qs2(dataset)
  if(!is.null(ts_data)) assign(dataset, ts_data, envir = .GlobalEnv)
  return(invisible(NULL))
}