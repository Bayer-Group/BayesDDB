#' Estimate propensity scores
#'
#' @keywords internal
#' @noRd
ps_predict <- function(model) {
  # ps_predict() estimates propensity scores
  # INPUT
  # model: fitted propensity score model
  # OUTPUT
  # Numeric vector of fitted probabilities
  prd <- tryCatch(
    stats::predict(model, type = "response"),
    error = function(e) {
      stop("Failed to compute propensity scores: ", e$message, call. = FALSE)
    }
  )

  if (!all(is.finite(prd))) {
    stop("Non-finite propensity scores produced by the PS model", call. = FALSE)
  }

  eps <- 1e-6
  if (any(prd < eps | prd > 1 - eps)) {
    warning(
      "Some propensity scores are numerically near 0 or 1; estimates may be unstable",
      call. = FALSE
    )
  }
  prd
}
