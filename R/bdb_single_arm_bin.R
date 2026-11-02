#' Bayesian Dynamic Borrowing for Single-Arm Trial with Binary Outcome
#'
#' Applies the Bayesian Dynamic Borrowing double-adjustment method for single-arm
#' trials with binary outcomes. The single-arm trial is augmented with external
#' data. No treatment variable is used.
#'
#' @param Y Numeric vector of binary outcomes (0s and 1s).
#' @param Z Numeric vector indicating data source: 1 for current trial data,
#'   0 for external data.
#' @param X Matrix or data frame of numeric baseline covariates (confounders)
#'   with n rows and p columns. Column names must be non-empty, unique and must
#'   not be one of \code{Y}, \code{Z}, \code{TRT}, \code{prd}, \code{stratum}
#'   or \code{gamma_discount}.
#' @param n_target_borrow Positive integer. Number of subjects intended to be
#'   borrowed from the external data source.
#' @param n_strata Positive integer. Number of strata of the propensity score
#'   distribution.
#' @param map_a Numeric scalar (finite, > 0). Tuning shape parameter for the arctangent elastic
#'   function. Larger values are more lenient towards outcome differences. See
#'   \code{\link{elastic_fun}}.
#' @param map_b Numeric scalar (finite, > 0). Tuning power parameter for the arctangent elastic
#'   function. See \code{\link{elastic_fun}}.
#' @param n_ppp Positive integer. Number of samples for estimating the posterior
#'   predictive probability (default 1000).
#' @param alpha_level Numeric in (0, 1). Significance level, e.g. 0.05 for 95\%
#'   intervals.
#' @param n_iter Positive integer. Number of samples from the posterior
#'   distribution.
#' @param alpha_shape Numeric scalar (finite, > 0). Shape1 hyperparameter of the Beta prior for
#'   the event probability.
#' @param beta_shape Numeric scalar (finite, > 0). Shape2 hyperparameter of the Beta prior for
#'   the event probability.
#' @param doubleadj Logical. If TRUE, both baseline and outcome adjustment. If
#'   FALSE, baseline adjustment only. Outcome adjustment is only applied with
#'   \code{risk_diff = FALSE}; see \code{risk_diff}.
#' @param risk_diff Logical. If FALSE, \code{theta} is the event rate of the
#'   single-arm trial pooled with the discounted external data. If TRUE,
#'   \code{theta} is the risk difference between the trial arm and the
#'   discounted external data (trial minus external). In this case the external
#'   data serve as the comparator, and their outcome is expected to differ from
#'   the treated trial arm by the treatment effect. The posterior predictive
#'   probability would then measure the treatment effect rather than prior-data
#'   conflict, so outcome adjustment is not applied: with
#'   \code{doubleadj = TRUE} a warning is issued and only baseline adjustment
#'   is used (\code{gPPP} is set to 1).
#' @param plots Logical. If TRUE, diagnostic plots are added to the output.
#' @param seed Optional non-negative whole number below \eqn{2^{31}}{2^31}
#'   used to seed the pseudo-random number generator (default NULL, no seed
#'   is set).
#'   The caller's random number generator state is restored on exit.
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
#'   \item{\code{ne_s1}, \code{ne_s0}}{Number of events in the current trial and
#'     in the external data by stratum.}
#'   \item{\code{r_s1}, \code{r_s0}}{\code{n_strata}-by-3 matrices with event
#'     rate and exact (Clopper-Pearson) confidence interval of the current trial
#'     and of the external data by stratum.}
#'   \item{\code{r_sat}, \code{r_eca}}{Overall event rate and confidence interval
#'     of the single-arm trial (SAT) and external control arm (ECA).}
#'   \item{\code{v}, \code{r}}{Overlapping coefficients and normalized
#'     overlapping coefficients by stratum.}
#'   \item{\code{A_borrow}, \code{A_eff}}{Effective number of subjects borrowed
#'     from the external data by stratum and overall.}
#'   \item{\code{PPP}, \code{gPPP}}{Posterior predictive probability and mapped
#'     posterior predictive probability by stratum. With
#'     \code{doubleadj = FALSE} or \code{risk_diff = TRUE}, \code{PPP} is still
#'     computed but \code{gPPP} is set to 1 (no outcome adjustment).}
#'   \item{\code{gamma_base}, \code{gamma_base_out}}{Baseline discounting factor
#'     and baseline and outcome discounting factor by stratum.}
#'   \item{\code{w_s}}{Weights for combining results from strata.}
#'   \item{\code{theta_m}}{\code{n_iter} pooled posterior samples.}
#'   \item{\code{post_mean_theta}, \code{post_var_theta}}{Posterior mean and
#'     variance.}
#'   \item{\code{CrI}, \code{Width}}{Credible interval and its length.}
#'   \item{\code{theta_s}}{\code{n_strata}-by-3 matrix of estimates and credible
#'     intervals by stratum.}
#'   \item{\code{theta_o}}{1-by-3 matrix of the overall estimate and credible
#'     interval.}
#' }
#' When \code{risk_diff = TRUE}, \code{theta} refers to the risk difference
#' rather than the pooled event rate. When \code{plots = TRUE}, the list also
#' contains the ggplot objects \code{estimate_plot}, \code{ps_group_plot},
#' \code{ps_group_stratum_plot}, \code{density_plot}, \code{overall_mean_plot}
#' and \code{stratum_mean_plot}.
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
#' @seealso \code{\link{bdb_single_arm_cont}}, \code{\link{bdb_two_arms_bin}},
#'   \code{\link{bdb_two_arms_cont}}, \code{\link{elastic_fun}},
#'   \code{\link{data_bin_sat}}
#'
#' @examples
#' \donttest{
#' # Single-arm trial with N1 = 200 subjects, borrowing information equivalent
#' # to 200 subjects from an external source of N0 = 2000 patients
#' out <- bdb_single_arm_bin(Y = data_bin_sat$Y,
#'                           Z = data_bin_sat$Z,
#'                           X = data_bin_sat[, paste0('X', 1:10)],
#'                           n_target_borrow = 200,
#'                           n_strata = 3,
#'                           map_a = 1,
#'                           map_b = 3,
#'                           n_ppp = 1e4,
#'                           alpha_level = 0.05,
#'                           n_iter = 1e5,
#'                           alpha_shape = 0.5,
#'                           beta_shape = 0.5,
#'                           doubleadj = TRUE,
#'                           risk_diff = FALSE,
#'                           plots = TRUE,
#'                           seed = 12345)
#' out$theta_o
#' }
#' @export
bdb_single_arm_bin <- function(
  Y,
  Z,
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
  risk_diff,
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
    TRT = NULL,
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
    risk_diff,
    plots,
    seed
  )

  ##############################################################################
  # Set pseudo-random generator seed
  local_seed(seed)

  ##############################################################################
  # Define data object
  data <- data.frame(Y, Z, X)
  ##############################################################################
  # Step 1: Model exposure mechanism
  ##############################################################################
  # Define exposure mechanism
  mylogit <- fit_ps_model(data, formula = Z ~ . - Y)

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

  # Indices current subjects
  is_current <- data_trim$Z == 1
  # Indices external subjects
  is_external <- data_trim$Z == 0

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

  # Current trial
  n_s1 <- stratum_counts(data_trim$stratum, is_current, n_strata)
  # n_1, overall number of subjects in the current trial
  n_1 <- sum(n_s1)
  # ne_s1, number of events in current trial by stratum
  ne_s1 <- stratum_events(data_trim$Y, data_trim$stratum, is_current, n_strata)
  # ne_1, overall number of events in the current trial
  ne_1 <- sum(ne_s1)

  # External Control
  # n_s0, number of subject in external control by stratum
  n_s0 <- stratum_counts(data_trim$stratum, is_external, n_strata)
  # n_0, overall number of subjects in the external control
  n_0 <- sum(n_s0)
  # ne_s0, number of event for external data by stratum
  ne_s0 <- stratum_events(data_trim$Y, data_trim$stratum, is_external, n_strata)
  # ne_0, overall number of events in the external in control
  ne_0 <- sum(ne_s0)

  # Weights for combining results from strata (through a weighted average of trial subjects)
  w_s <- n_s1 / sum(n_s1)

  # There must be at least two subjects on each stratum otherwise calculations cannot be computed
  if (any(n_s0 < 2) || any(n_s1 < 2)) {
    stop(
      "There must be at least two subjects on each stratum in both the current trial ",
      "and the external data. Consider reducing n_strata or revising the propensity ",
      "score model.",
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
    ne_ref = ne_s1,
    n_ref = n_s1,
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

  # With risk_diff = TRUE the external data are the comparator, so the PPP
  # compares the treated trial arm with the external controls and reflects the
  # treatment effect itself. Outcome adjustment is therefore not applied.
  if (risk_diff && doubleadj) {
    gPPP <- rep(1, n_strata)
    warning(
      "Double adjustment is not valid for the risk difference estimator. ",
      "Outcome adjustment is not considered; only baseline adjustment is applied.",
      call. = FALSE
    )
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

  # Define baseline and outcome discount parameter
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

  if (risk_diff) {
    # Initialize matrices for later use
    theta <- matrix(0, nrow = n_strata, ncol = n_iter)
    theta_s <- matrix(0, nrow = n_strata, ncol = 3)
    colnames(theta_s) <- c("mean", "lower", "upper")

    for (i in 1:n_strata) {
      # # Treatment arm
      # Sampling from beta-binomial distribution with parameters shape parameters alpha_shape and beta_shape
      thetaT_i <- stats::rbeta(
        n_iter,
        alpha_shape + ne_s1[i],
        beta_shape + n_s1[i] - ne_s1[i]
      )
      # Controls
      # Sampling from beta-binomial distribution with parameters shape parameters alpha_shape and beta_shape
      thetaC_i <- stats::rbeta(
        n_iter,
        alpha_shape + ne_s0[i] * gamma_base_out[i],
        beta_shape + (n_s0[i] - ne_s0[i]) * gamma_base_out[i]
      )
      # Treatment effect (Treatment - Controls)
      theta[i, ] <- thetaT_i - thetaC_i
      # Mean rate by stratum
      theta_s[i, ] <- draws_summary(theta[i, ], alpha_level)
      # Discount parameter
      idx_CONT <- is_external & data_trim$stratum == i
      data_trim$gamma_discount[idx_CONT] <- gamma_base_out[i]
    }
  } else {
    # Initialize matrices for later use
    theta <- matrix(0, nrow = n_strata, ncol = n_iter)
    theta_s <- matrix(0, nrow = n_strata, ncol = 3)
    colnames(theta_s) <- c("mean", "lower", "upper")
    for (i in 1:n_strata) {
      # Single-arm trial with discounted external data
      theta[i, ] <- stats::rbeta(
        n_iter,
        alpha_shape + ne_s1[i] + ne_s0[i] * gamma_base_out[i],
        beta_shape +
          (n_s1[i] - ne_s1[i]) +
          (n_s0[i] - ne_s0[i]) * gamma_base_out[i]
      )
      # Posterior mean by stratum and 95% CrI
      theta_s[i, ] <- draws_summary(theta[i, ], alpha_level)
      # Discount parameter
      idx_CONT <- is_external & data_trim$stratum == i
      data_trim$gamma_discount[idx_CONT] <- gamma_base_out[i]
    }
  }

  pool_results <- pool_strata(theta, w_s, alpha_level)

  ##############################################################################
  ##############################################################################
  # Current treatment event rate
  ##############################################################################
  ##############################################################################
  # Event rate by stratum and CI (Clopper-Pearson)
  r_s1 <- binom_ci_exact(ne_s1, n_s1, alpha_level)

  # Overall event rate and CI from SAT (Clopper-Pearson)
  r_sat <- binom_ci_exact(ne_1, n_1, alpha_level)

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
  # Function output
  output <- list(
    data_trim = data_trim,
    N0_trim = n_external,
    ps_model = mylogit,
    ps_breaks = q,
    n_s1 = n_s1,
    n_s0 = n_s0,
    ne_s1 = ne_s1,
    ne_s0 = ne_s0,
    r_s1 = r_s1,
    r_s0 = r_s0,
    r_sat = r_sat,
    r_eca = r_eca,
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
      title = if (risk_diff) {
        "Posterior Risk Difference by Stratum and Overall"
      } else {
        "Posterior SAT Pooled Event Rate by Stratum and Overall"
      },
      digits = 5,
      hline = risk_diff
    )

    # Plot - Propensity Score Distributions by Group
    ps_group_plot <- plot_ps_by_group(data_trim$prd, data_trim$Z)

    # Plot - Propensity Score Distributions by Stratum
    ps_group_stratum_plot <- plot_ps_by_stratum(
      data_trim$prd,
      data_trim$stratum,
      data_trim$Z
    )

    # Plot - Posterior Density
    density_plot <- plot_posterior_density(pool_results$theta_m)

    # Plot - Incidence proportion for ECA and SAT
    overall_mean_plot <- plot_group_overall(
      external = r_eca,
      internal = r_sat,
      internal_label = "SAT",
      alpha_level = alpha_level,
      value_label = "Event rate",
      title = "Overall event rate for ECA and SAT"
    )

    # Plot - Forest plot of ICA and ECA event rates by stratum
    stratum_mean_plot <- plot_group_by_stratum(
      external = r_s0,
      internal = r_s1,
      internal_label = "SAT",
      n_strata = n_strata,
      alpha_level = alpha_level,
      value_label = "Event rate",
      title = "Event rate for ECA and SAT by stratum"
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
  return(output)
  # function ends
}
