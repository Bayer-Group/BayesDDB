#' Posterior predictive probability for a binary outcome
#'
#' @keywords internal
#' @noRd
ppp_binary <- function(
  ne_s0,
  n_s0,
  ne_ref,
  n_ref,
  alpha_shape,
  beta_shape,
  n_ppp
) {
  # Posterior predictive probability for a binary outcome
  #
  # Draws event probabilities from the beta posterior of the external data, then
  # simulates the reference arm under those probabilities and compares against
  # its observed event rate. The reported value is the mid-p average of the
  # strict and non-strict comparisons.
  #
  # INPUT
  # ne_s0: External events per stratum.
  # n_s0: External subjects per stratum.
  # ne_ref: Reference-arm events per stratum.
  # n_ref: Reference-arm subjects per stratum.
  # alpha_shape: Beta prior hyper-parameter.
  # beta_shape: Beta prior hyper-parameter.
  # n_ppp: Number of simulated samples.
  # OUTPUT
  # Double vector of probabilities, one per stratum.
  ################################################################################
  n_strata <- length(n_s0)

  # Observed mean in the reference arm by stratum
  e.m <- ne_ref / n_ref

  # beta-binomial posterior for external data by stratum
  prb <- matrix(0, nrow = n_strata, ncol = n_ppp)
  for (i in seq_len(n_strata)) {
    prb[i, ] <- stats::rbeta(
      n_ppp,
      shape1 = alpha_shape + ne_s0[i],
      shape2 = beta_shape + n_s0[i] - ne_s0[i]
    )
  }

  # Posterior predictive probability (sampling from current study with
  # probabilities from external data and comparing against the current study)
  PPP <- rep(0, n_strata)
  for (i in seq_len(n_strata)) {
    p_one_strict <- mean(
      stats::rbinom(n_ppp, n_ref[i], prb[i, ]) / n_ref[i] < e.m[i]
    )
    p_one_non_strict <- mean(
      stats::rbinom(n_ppp, n_ref[i], prb[i, ]) / n_ref[i] <= e.m[i]
    )
    # use mid-p value for bias correction
    PPP[i] <- (p_one_strict + p_one_non_strict) / 2
  }
  PPP
}
