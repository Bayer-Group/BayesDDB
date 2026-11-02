#' Bayesian Dynamic Borrowing for Two-Arm Trial with Binary Outcome
#'
#' Applies the Bayesian Dynamic Borrowing double-adjustment method when the control arm
#' of a randomized trial is augmented with external data (binary outcome).
#'
#' @param Y Numeric vector of binary outcomes (0s and 1s).
#' @param Z Numeric vector indicating data source: 1 for current trial data, 0 for external data.
#' @param TRT Numeric vector indicating treatment: 1 for active treatment, 0 for
#'   control. External subjects must be controls.
#' @param doubleadj Logical. If TRUE, both baseline and outcome adjustment. If
#'   FALSE, baseline adjustment only.
#' @inheritParams bdb_single_arm_bin
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{data_trim}}{Trimmed data set including the estimated propensity
#'     score (\code{prd}), \code{stratum} and discount value
#'     (\code{gamma_discount}) per subject.}
#'   \item{\code{N0_trim}}{Overall number of external subjects after trimming.}
#'   \item{\code{ps_model}}{Fitted propensity score model (\code{glm} object).}
#'   \item{\code{ps_breaks}}{Stratum cut-points for the propensity scores.}
#'   \item{\code{n_s1}, \code{n_s0}}{Number of current and external subjects by
#'     stratum.}
#'   \item{\code{n_s11}, \code{n_s10}}{Number of current subjects in the
#'     treatment arm and in the internal control arm by stratum.}
#'   \item{\code{ne_s11}, \code{ne_s10}, \code{ne_s0}}{Number of events in the
#'     treatment arm, internal control arm and external data by stratum.}
#'   \item{\code{r_s11}, \code{r_s10}, \code{r_s0}}{\code{n_strata}-by-3 matrices
#'     with event rate and exact (Clopper-Pearson) confidence interval of the
#'     treatment arm, internal control arm and external data by stratum.}
#'   \item{\code{r_ica}, \code{r_eca}}{Overall event rate and confidence interval
#'     of the internal control arm (ICA) and external control arm (ECA).}
#'   \item{\code{theta_n}}{Treatment effect from the current trial alone (without
#'     using external data) with credible interval.}
#'   \item{\code{v}, \code{r}}{Overlapping coefficients and normalized
#'     overlapping coefficients by stratum.}
#'   \item{\code{A_borrow}, \code{A_eff}}{Effective number of subjects borrowed
#'     from the external data by stratum and overall.}
#'   \item{\code{PPP}, \code{gPPP}}{Posterior predictive probability and mapped
#'     posterior predictive probability by stratum. With
#'     \code{doubleadj = FALSE}, \code{PPP} is still computed but \code{gPPP}
#'     is set to 1 (no outcome adjustment).}
#'   \item{\code{gamma_base}, \code{gamma_base_out}}{Baseline discounting factor
#'     and baseline and outcome discounting factor by stratum.}
#'   \item{\code{w_s}}{Weights for combining results from strata.}
#'   \item{\code{theta_m}}{\code{n_iter} pooled posterior samples of the
#'     treatment effect.}
#'   \item{\code{post_mean_theta}, \code{post_var_theta}}{Posterior mean and
#'     variance of the treatment effect (difference of event rates between
#'     treatment and control arms).}
#'   \item{\code{CrI}, \code{Width}}{Credible interval and its length.}
#'   \item{\code{theta_s}}{\code{n_strata}-by-3 matrix of treatment effect
#'     estimates and credible intervals by stratum.}
#'   \item{\code{theta_o}}{1-by-3 matrix of the overall treatment effect and
#'     credible interval.}
#' }
#' When \code{plots = TRUE}, the list also contains the ggplot objects
#' \code{estimate_plot}, \code{ps_group_plot}, \code{ps_group_stratum_plot},
#' \code{density_plot}, \code{overall_mean_plot} and \code{stratum_mean_plot}.
#'
#' @references
#' Farjat, A., Ji, Y., Kaiser, A., Lu, C., Pap, A., Potts, J., Wang, M.
#' A Comprehensive Bayesian Double-Adjustment Approach to Dynamic Borrowing of
#' External Data. Manuscript in progress.
#'
#' Wang, C., Li, H., Chen, W.C., Lu, N., Tiwari, R., Xu, Y. and Yue, L.Q., 2019.
#' Propensity score-integrated power prior approach for incorporating real-world
#' evidence in single-arm clinical studies.
#' Journal of biopharmaceutical statistics, 29(5), pp.731-748.
#'
#' Jiang, L., Nie, L. and Yuan, Y., 2023. Elastic priors to dynamically borrow
#' information from historical data in clinical trials. Biometrics, 79(1),
#' pp.49-60. \doi{10.1111/biom.13551}
#'
#' @seealso \code{\link{bdb_single_arm_bin}}, \code{\link{bdb_single_arm_cont}},
#'   \code{\link{bdb_two_arms_cont}}, \code{\link{elastic_fun}},
#'   \code{\link{data_bin}}
#'
#' @examples
#' \donttest{
#' # Trial with N1 = 300 subjects randomized 2:1 (treatment vs control),
#' # borrowing information equivalent to 100 subjects from an external source
#' # of N0 = 3000 patients to achieve a 1:1 ratio of treatment to controls
#' out <- bdb_two_arms_bin(Y = data_bin$Y,
#'                         Z = data_bin$Z,
#'                         TRT = data_bin$TRT,
#'                         X = data_bin[, paste0('X', 1:10)],
#'                         n_target_borrow = 100,
#'                         n_strata = 3,
#'                         map_a = 1,
#'                         map_b = 1,
#'                         n_ppp = 1e4,
#'                         alpha_level = 0.05,
#'                         n_iter = 1e5,
#'                         alpha_shape = 0.5,
#'                         beta_shape = 0.5,
#'                         doubleadj = TRUE,
#'                         plots = TRUE,
#'                         seed = 12345)
#' out$theta_o
#' }
#' @export
bdb_two_arms_bin <- function(
  Y,
  Z,
  TRT,
  X,
  n_target_borrow,
  n_strata,
  map_a,
  map_b,
  n_ppp = 1e3,
  alpha_level,
  n_iter,
  alpha_shape,
  beta_shape,
  doubleadj,
  plots,
  seed = NULL
) {
  ##############################################################################
  # Function begins
  # Check input arguments
  input_checks_bin(
    Y,
    Z,
    X,
    TRT,
    n_target_borrow,
    n_strata,
    map_a,
    map_b,
    n_ppp,
    alpha_level,
    n_iter,
    alpha_shape,
    beta_shape,
    doubleadj,
    risk_diff = NULL,
    plots,
    seed
  )

  ##############################################################################
  # Set pseudo-random generator seed
  local_seed(seed)

  ##############################################################################
  # Define data object
  data <- data.frame(Y, Z, TRT, X)
  ##############################################################################
  # Step 1: Model exposure mechanism
  ##############################################################################
  # Define exposure mechanism
  mylogit <- fit_ps_model(data, formula = Z ~ . - Y - TRT)

  ##############################################################################
  ##############################################################################
  # Step 2: Propensity Score Calculation
  ##############################################################################
  ##############################################################################
  # Propensity scores
  prd <- ps_predict(mylogit)

  # Add propensity scores to data object
  data <- data.frame(data, prd)

  # Propensity scores of current patients
  prd_current <- prd[Z == 1]

  ##############################################################################
  ##############################################################################
  # Step 3: Trimming
  ##############################################################################
  ##############################################################################
  # Define data_trim and restrict propensity scores
  # to the range of values of subjects in the current trial
  data_trim <- trimming(data, prd, prd_current, sort_by_z = TRUE)

  is_current <- data_trim$Z == 1
  is_external <- data_trim$Z == 0
  is_treated <- data_trim$Z == 1 & data_trim$TRT == 1
  is_int_control <- data_trim$Z == 1 & data_trim$TRT == 0

  # Trimmed Propensity scores
  # Propensity score of patients in the trial
  prd1 <- data_trim$prd[is_current]

  # Overall number of subjects after trimming
  n <- nrow(data_trim)

  # Number of external subjects after trimming
  n_external <- sum(is_external)

  # Checks after trimming
  if (n == 0) {
    stop("No observations remain after trimming")
  }
  if (n_external == 0) {
    stop("No external observations remain after trimming")
  }

  ##############################################################################
  ##############################################################################
  # Step 4: Stratification
  ##############################################################################
  ##############################################################################
  # Strata and PS breaks
  strata <- stratify_ps(data_trim$prd, prd1, n_strata)
  data_trim$stratum <- strata$stratum
  q <- strata$breaks

  # Current trial (both arms)
  # n_s1, number subjects in the current trial by stratum
  n_s1 <- stratum_counts(data_trim$stratum, is_current, n_strata)

  # Treatment arm
  # n_s11, number subjects in the treatment arm by stratum
  n_s11 <- stratum_counts(data_trim$stratum, is_treated, n_strata)
  # n_11, overall number of subjects in the treatment arm
  n_11 <- sum(n_s11)
  # ne_s11, number of events in the treatment arm by stratum
  ne_s11 <- stratum_events(data_trim$Y, data_trim$stratum, is_treated, n_strata)
  # ne_11, overall number of events in the current trial
  ne_11 <- sum(ne_s11)

  # Internal Control
  # n_s10, number subjects in the internal control by stratum
  n_s10 <- stratum_counts(data_trim$stratum, is_int_control, n_strata)
  # n_10, overall number of subjects in the internal control arm
  n_10 <- sum(n_s10)
  # ne_s10, number of events in the internal control arm by stratum
  ne_s10 <- stratum_events(
    data_trim$Y,
    data_trim$stratum,
    is_int_control,
    n_strata
  )
  # ne_10, overall number of events in the internal control arm
  ne_10 <- sum(ne_s10)

  # External Control
  # n_s0, number of subject in external control arm by stratum
  n_s0 <- stratum_counts(data_trim$stratum, is_external, n_strata)
  # n_0, overall number of subjects in the external control arm
  n_0 <- sum(n_s0)
  # ne_s0, number of event in the external control arm by stratum
  ne_s0 <- stratum_events(data_trim$Y, data_trim$stratum, is_external, n_strata)
  # ne_0, overall number of events in the external control arm
  ne_0 <- sum(ne_s0)

  # Weights for combining results from different strata (through a weighted average of trial subjects)
  w_s <- n_s1 / sum(n_s1)

  # There must be at least two subjects on each cell otherwise calculations
  # cannot be computed
  if (any(n_s11 < 2) || any(n_s10 < 2) || any(n_s0 < 2)) {
    stop(
      "There must be at least two subjects on each stratum in the treatment arm, the ",
      "internal control arm and the external data. Consider reducing n_strata or ",
      "revising the propensity score model.",
      call. = FALSE
    )
  }

  ##############################################################################
  ##############################################################################
  # Step 5: Overlapping coefficient
  ##############################################################################
  ##############################################################################
  # Calculate vector of overlapping coefficients
  v <- overlap_by_stratum(
    data_trim$prd,
    data_trim$Z,
    data_trim$stratum,
    n_strata
  )

  ##############################################################################
  ##############################################################################
  # Step 6: Weighting
  ##############################################################################
  ##############################################################################
  # Initiate vector of normalized overlapping coefficients
  r <- v / sum(v)

  ##############################################################################
  ##############################################################################
  # Step 7: Calibration
  ##############################################################################
  ##############################################################################
  PPP <- ppp_binary(
    ne_s0 = ne_s0,
    n_s0 = n_s0,
    ne_ref = ne_s10,
    n_ref = n_s10,
    alpha_shape = alpha_shape,
    beta_shape = beta_shape,
    n_ppp = n_ppp
  )

  ##############################################################################
  ##############################################################################
  # Step 8: Mapping
  ##############################################################################
  ##############################################################################
  # Mapped predictive posterior probability
  gPPP <- elastic_fun(PPP, map_a, map_b)

  # In case of baseline adjustment only
  if (!doubleadj) {
    gPPP <- rep(1, n_strata)
  }

  ##############################################################################
  ##############################################################################
  # Step 9: Discounting
  ##############################################################################
  ##############################################################################
  # Define baseline discounting factor
  gamma_base <- rep(0, n_strata)
  for (i in 1:n_strata) {
    gamma_base[i] <- min(1, n_target_borrow * r[i] / n_s0[i])
  }

  # Define baseline and outcome discounting factor
  gamma_base_out <- rep(0, n_strata)
  for (i in 1:n_strata) {
    gamma_base_out[i] <- min(1, n_target_borrow * r[i] * gPPP[i] / n_s0[i])
  }

  if (n_target_borrow > sum(n_s0)) {
    warning(
      "Requested borrowing n_target_borrow exceeds number of trimmed external subjects"
    )
  }

  # Number of subjects actually borrowed from external source by stratum
  A_borrow <- n_s0 * gamma_base_out

  # Effective sample size
  A_eff <- sum(A_borrow)

  # Initialize vector of discount values (current trial: 1, external control: gamma_base_out )
  data_trim$gamma_discount <- 1

  ##############################################################################
  ##############################################################################
  # Steps 10-11: Analysis/Summary
  ##############################################################################
  ##############################################################################

  ##############################################################################
  ##############################################################################
  # Both baseline and outcome adjustment
  ##############################################################################
  ##############################################################################
  # Initialize matrices for later use
  theta <- matrix(0, nrow = n_strata, ncol = n_iter)
  theta_s <- matrix(0, nrow = n_strata, ncol = 3)
  colnames(theta_s) <- c("mean", "lower", "upper")
  for (i in 1:n_strata) {
    # # Treatment arm
    # Sampling from beta-binomial distribution with parameters shape parameters alpha_shape and beta_shape
    thetaT_i <- stats::rbeta(
      n_iter,
      alpha_shape + ne_s11[i],
      beta_shape + n_s11[i] - ne_s11[i]
    )

    # Controls
    # Sampling from beta-binomial distribution with parameters shape parameters alpha_shape and beta_shape
    thetaC_i <- stats::rbeta(
      n_iter,
      alpha_shape + ne_s10[i] + ne_s0[i] * gamma_base_out[i],
      beta_shape +
        (n_s10[i] - ne_s10[i]) +
        (n_s0[i] - ne_s0[i]) * gamma_base_out[i]
    )

    # Treatment effect (Treatment - Controls)
    theta[i, ] <- thetaT_i - thetaC_i
    # Mean rate by stratum
    theta_s[i, ] <- draws_summary(theta[i, ], alpha_level)
    # Discount parameter
    idx_CONT <- data_trim$TRT == 0 & data_trim$stratum == i
    z_s <- data_trim$Z[idx_CONT]
    gamma_s <- ifelse(z_s == 1, 1, gamma_base_out[i])
    data_trim$gamma_discount[idx_CONT] <- gamma_s
  }

  pool_results <- pool_strata(theta, w_s, alpha_level)

  ##############################################################################
  ##############################################################################
  # Current treatment event rate
  ##############################################################################
  ##############################################################################
  # Event rate by stratum and CI (Clopper-Pearson)
  r_s11 <- binom_ci_exact(ne_s11, n_s11, alpha_level)

  ##############################################################################
  ##############################################################################
  # ICA event rate
  ##############################################################################
  ##############################################################################
  # Event rate by stratum and CI (Clopper-Pearson)
  r_s10 <- binom_ci_exact(ne_s10, n_s10, alpha_level)

  # Overall event rate and CI from ICA (Clopper-Pearson)
  r_ica <- binom_ci_exact(ne_10, n_10, alpha_level)

  ##############################################################################
  ##############################################################################
  # ECA event rate
  ##############################################################################
  ##############################################################################
  # Event rate by stratum and CI (Clopper-Pearson)
  r_s0 <- binom_ci_exact(ne_s0, n_s0, alpha_level)

  # Overall event rate and CI from ECA (Clopper-Pearson)
  r_eca <- binom_ci_exact(ne_0, n_0, alpha_level)

  ##############################################################################
  ##############################################################################
  # Treatment effect from the current trial (without using external data)
  # Difference of rates (treatment - control)
  diff_n <- stats::rbeta(
    n_iter,
    alpha_shape + ne_11,
    beta_shape + n_11 - ne_11
  ) -
    stats::rbeta(n_iter, alpha_shape + ne_10, beta_shape + n_10 - ne_10)
  # Posterior mean and credible interval
  post.mean_n <- mean(diff_n)
  CrI_n <- as.numeric(stats::quantile(
    diff_n,
    c(alpha_level / 2, 1 - alpha_level / 2)
  ))
  theta_n <- c(post.mean_n, CrI_n)
  names(theta_n) <- c("mean", "lower", "upper")

  ##############################################################################
  # Function output
  output <- list(
    data_trim = data_trim,
    N0_trim = n_external,
    ps_model = mylogit,
    ps_breaks = q,
    n_s1 = n_s1,
    n_s0 = n_s0,
    n_s11 = n_s11,
    n_s10 = n_s10,
    ne_s11 = ne_s11,
    ne_s10 = ne_s10,
    ne_s0 = ne_s0,
    r_s11 = r_s11,
    r_s10 = r_s10,
    r_s0 = r_s0,
    r_ica = r_ica,
    r_eca = r_eca,
    theta_n = theta_n,
    v = v,
    r = r,
    A_borrow = A_borrow,
    A_eff = A_eff,
    PPP = PPP,
    gPPP = gPPP,
    gamma_base = gamma_base,
    gamma_base_out = gamma_base_out,
    w_s = w_s,
    theta_m = pool_results$theta_m,
    post_mean_theta = pool_results$post_mean_theta,
    post_var_theta = pool_results$post_var_theta,
    CrI = pool_results$CrI,
    Width = pool_results$Width,
    theta_s = theta_s,
    theta_o = pool_results$theta_o
  )

  ##############################################################################
  if (plots) {
    # Plot - Forest plot of treatment effect by stratum and overall
    estimate_plot <- plot_forest_strata(
      theta_s = theta_s,
      theta_o = pool_results$theta_o,
      n_strata = n_strata,
      alpha_level = alpha_level,
      title = "Treatment Effect by Stratum and Overall",
      digits = 5,
      hline = TRUE
    )

    # Plot - Propensity Score Distributions by Group
    ps_group_plot <- plot_ps_by_group(data_trim$prd, data_trim$Z)

    # Plot - Propensity Score Distributions by Stratum
    ps_group_stratum_plot <- plot_ps_by_stratum(
      data_trim$prd,
      data_trim$stratum,
      data_trim$Z
    )

    # Plot - Posterior density
    density_plot <- plot_posterior_density(pool_results$theta_m)

    # Plot - Incidence proportion for ECA and ICA
    overall_mean_plot <- plot_group_overall(
      external = r_eca,
      internal = r_ica,
      internal_label = "ICA",
      alpha_level = alpha_level,
      value_label = "Event rate",
      title = "Overall event rate for ECA and ICA"
    )

    # Plot - Forest plot of ICA and ECA event rates by stratum
    stratum_mean_plot <- plot_group_by_stratum(
      external = r_s0,
      internal = r_s10,
      internal_label = "ICA",
      n_strata = n_strata,
      alpha_level = alpha_level,
      value_label = "Event rate",
      title = "Event rate for ECA and ICA by stratum"
    )

    # Function output
    output <- c(
      output,
      list(
        estimate_plot = estimate_plot,
        ps_group_plot = ps_group_plot,
        ps_group_stratum_plot = ps_group_stratum_plot,
        density_plot = density_plot,
        overall_mean_plot = overall_mean_plot,
        stratum_mean_plot = stratum_mean_plot
      )
    )
  }
  ##############################################################################
  ##############################################################################
  return(output)
  # function ends
}
