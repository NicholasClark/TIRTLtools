#' @usage TIRTL_process(...)
#' @keywords internal
TIRTL_process <- function(...) {
  lifecycle::deprecate_warn(when = "0.2.1", what="TIRTL_process()",
                            details = "This function has been deprecated because its functionality has been added to the load_tirtlseq() function"
  )
  process_TIRTLseq(...)
}
