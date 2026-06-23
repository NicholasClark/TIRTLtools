#' Load example TIRTL-seq data
#'
#' @description
#' `r lifecycle::badge('experimental')`
#' This is a helper function to load example TIRTL-seq data.
#' The function currently loads one of two available TIRTL-seq datasets. It saves the data
#' to the global environment in an object with the same name as the dataset.
#' For example, if \code{dataset} is "SJTRC_minimal", it will load the data to an object
#' named "SJTRC_minimal" in the global environment. If the dataset is already loaded, it will not be re-loaded.
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
#' The function returns NULL invisibly, but adds a
#'
#' @family data_loading
#'
#' @examples
#' load_example_data("SJTRC_minimal")
#' load_example_data("SJTRC_longitudinal")
#'
#'

load_example_data = function(dataset = c("SJTRC_minimal", "SJTRC_longitudinal"), verbose = TRUE) {
  dataset = dataset[1]
  assert_choice(dataset, choices = c("SJTRC_minimal", "SJTRC_longitudinal"))

  if(dataset == "SJTRC_minimal") {
    if (exists("SJTRC_minimal", envir = .GlobalEnv, inherits = FALSE)) {
      if(verbose) message("Example data already loaded into object: 'SJTRC_minimal'")
      return(invisible(NULL))
    }
    ts_data = load_sjtrc_minimal()
    assign("SJTRC_minimal", ts_data, envir = .GlobalEnv)

    if(verbose) message("Example data loaded into object: 'SJTRC_minimal'")
  } else if(dataset == "SJTRC_longitudinal") {
    if (exists("SJTRC_longitudinal", envir = .GlobalEnv, inherits = FALSE)) {
      if(verbose) message("Example data already loaded into object: 'SJTRC_longitudinal'")
      return(invisible(NULL))
    }
    ts_data = load_sjtrc_longitudinal()
    assign("SJTRC_longitudinal", ts_data, envir = .GlobalEnv)
    if(verbose) message("Example data loaded into object: 'SJTRC_longitudinal'")
  }
  return(invisible(NULL))
}

load_sjtrc_minimal = function() {
  folder = system.file("extdata/SJTRC_TIRTLseq_minimal",
                       package = "TIRTLtools")
  ts_data = load_tirtlseq(folder,
                          meta_columns = c("marker", "timepoint", "version"), sep = "_")
  return(ts_data)
}

load_sjtrc_longitudinal = function() {
  folder = system.file("extdata/SJTRC_TIRTL_seq_longitudinal", package = "TIRTLtools")
  ts_data = load_tirtlseq(folder, meta_columns = c("marker", "timepoint", "version"), sep = "_", verbose = FALSE)
  return(ts_data)
}
