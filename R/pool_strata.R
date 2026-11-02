#' Pool per-stratum posterior draws
#'
#' @keywords internal
#' @noRd
pool_strata <- function(theta, w_s, alpha_level) {
  # Pool per-stratum posterior draws into an overall posterior
  # Takes the stratum-weighted mean of the draws at each iteration, then
  # summarizes the resulting sample.
  #
  # INPUT
  # theta: n_strata-by-n_iter matrix of posterior draws.
  # w_s: stratum weights.
  # alpha_level: Significance level.
  #
  # OUTPUT
  # theta_m: weighted mean across strata
  # post_mean_theta: posterior mean
  # post_var_theta: posterior variance
  # CrI: credible interval
  # Width: length of credible interval
  # theta_o: overall posterior mean with credible interval
  ##########################################################
  # Calculate weighted mean across strata (posterior distribution).
  theta_m <- as.vector(base::crossprod(w_s / sum(w_s), theta))
  post_mean_theta <- mean(theta_m)
  post_var_theta <- stats::var(theta_m)
  # Credible interval
  CrI <- as.numeric(stats::quantile(
    theta_m,
    c(alpha_level / 2, 1 - alpha_level / 2)
  ))
  # Length of credible interval
  Width <- CrI[2] - CrI[1]
  # Overall posterior mean with credible interval
  theta_o <- cbind(mean = post_mean_theta, lower = CrI[1], upper = CrI[2])

  list(
    theta_m = theta_m,
    post_mean_theta = post_mean_theta,
    post_var_theta = post_var_theta,
    CrI = CrI,
    Width = Width,
    theta_o = theta_o
  )
}
