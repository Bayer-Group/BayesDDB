#' Sample mean and t-based confidence interval
#'
#' @keywords internal
#' @noRd
# Function for estimating mean and confidence interval
mean_CI <- function(x, alpha_level) {
  # mean_CI() estimates mean and 100(1-alpha_level)% confidence interval
  # INPUT
  # x: array
  # alpha_level: significance level
  # OUTPUT
  # 1-by-3 matrix with mean, lower and upper 100(1-alpha_level)% confidence interval
  n <- length(x)
  if (n < 2) {
    stop("Need at least two observations to compute CI")
  }
  mu <- base::mean(x)
  sigma <- stats::sd(x)
  ta <- stats::qt(1 - alpha_level / 2, df = n - 1)
  output <- cbind(
    mean = mu,
    lower = mu - ta * sigma / sqrt(n),
    upper = mu + ta * sigma / sqrt(n)
  )
  return(output)
}
