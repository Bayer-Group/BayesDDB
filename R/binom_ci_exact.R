#' Event rate with exact (Clopper-Pearson) confidence interval
#'
#' @keywords internal
#' @noRd
binom_ci_exact <- function(x, n, alpha_level) {
  # binom_ci_exact() estimates event rate and exact confidence interval (Clopper-Pearson)
  # INPUT
  # x: number of events
  # n: total sample
  # alpha_level: significance level
  # OUTPUT
  # 1-by-3 matrix with mean, lower and upper 100(1-alpha_level)% confidence interval
  if (any(n <= 0)) {
    stop("binom_ci_exact requires n > 0")
  }
  output <- cbind(
    mean = x / n,
    lower = ifelse(x == 0, 0, stats::qbeta(alpha_level / 2, x, n - x + 1)),
    upper = ifelse(x == n, 1, stats::qbeta(1 - alpha_level / 2, x + 1, n - x))
  )
  return(output)
}
