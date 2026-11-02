#' Check whether a scalar is a whole number
#'
#' @keywords internal
#' @noRd
is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  # Check wholenumerness accounting for floating point precision
  # INPUT
  # x: scalar
  # OUTPUT
  # single logical value: TRUE (number is whole), FALSE (number is not whole)
  is.numeric(x) && length(x) == 1L && is.finite(x) && abs(x - round(x)) < tol
}
