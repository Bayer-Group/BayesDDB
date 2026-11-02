#' Posterior mean and credible interval of posterior draws
#'
#' @keywords internal
#' @noRd
draws_summary <- function(draws, alpha_level) {
  # Posterior mean and 100(1-alpha_level)% credible interval
  # INPUT
  # draws: Posterior draws
  # alpha_level: Significance level.
  # OUTPUT
  # 1-by-3 matrix with mean, lower and upper 100(1-alpha_level)% credible interval
  output <- cbind(
    mean = mean(draws),
    lower = as.numeric(stats::quantile(draws, alpha_level / 2)),
    upper = as.numeric(stats::quantile(draws, 1 - alpha_level / 2))
  )
  return(output)
}
