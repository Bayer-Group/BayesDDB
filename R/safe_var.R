#' Variance bounded away from zero
#'
#' @keywords internal
#' @noRd
# Estimates variance and prevents zero as outcome
# INPUT
# x: array
# OUTPUT
# variance
safe_var <- function(x, eps = 1e-8) {
  v <- stats::var(x)
  if (is.na(v) || v < eps) eps else v
}
