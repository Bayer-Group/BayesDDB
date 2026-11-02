#' Bayesian Dynamic Borrowing for Two-Arm Trial with Continuous Outcome
#'
#' Applies the Bayesian Dynamic Borrowing double-adjustment method when the control arm
#' of a randomized trial is augmented with external data (continuous outcome). Requires JAGS.
#'
#' @param Y Numeric vector of continuous outcomes.
#' @param Z Numeric vector indicating data source: 1 for current trial data, 0 for external data.
#' @param TRT Numeric vector indicating treatment: 1 for active treatment, 0 for
#'   control. External subjects must be controls.
#' @param doubleadj Logical. If TRUE, both baseline and outcome adjustment. If
#'   FALSE, baseline adjustment only.
#' @inheritParams bdb_single_arm_cont
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{data_trim}}{Trimmed data set including the estimated propensity
#'     score (\code{prd}), \code{stratum} and discount value
#'     (\code{gamma_discount}) per subject.}
#'   \item{\code{N0_trim}}{Overall number of external subjects after trimming.}
#'   \item{\code{ps_model}}{Fitted propensity score model (\code{glm} object).}
#'   \item{\code{ps_breaks}}{Stratum cut-points for the propensity scores.}
#'   \item{\code{n_s0}, \code{n_s1}}{Number of external and current subjects by
#'     stratum.}
#'   \item{\code{n_s11}, \code{n_s10}}{Number of current subjects in the
#'     treatment arm and in the internal control arm by stratum.}
#'   \item{\code{r_s0}, \code{r_s10}}{\code{n_strata}-by-3 matrices with sample
#'     mean and t-based confidence interval of the external control arm (ECA)
#'     and internal control arm (ICA) by stratum.}
#'   \item{\code{r_ica}, \code{r_eca}}{Overall sample mean and confidence
#'     interval of the ICA and ECA.}
#'   \item{\code{A_borrow}, \code{A_eff}}{Effective number of subjects borrowed
#'     from the external data by stratum and overall.}
#'   \item{\code{v}, \code{r}}{Overlapping coefficients and normalized
#'     overlapping coefficients by stratum.}
#'   \item{\code{gamma_base}, \code{gamma_base_out}}{Baseline discounting factor
#'     and baseline and outcome discounting factor by stratum.}
#'   \item{\code{PPP}, \code{gPPP}}{Posterior predictive probability and mapped
#'     posterior predictive probability by stratum. With
#'     \code{doubleadj = FALSE}, \code{PPP} is still computed but \code{gPPP}
#'     is set to 1 (no outcome adjustment).}
#'   \item{\code{w_s}}{Weights for combining results from strata.}
#'   \item{\code{theta_m}}{\code{n_iter} pooled posterior samples of the
#'     treatment effect.}
#'   \item{\code{post_mean_theta}, \code{post_var_theta}}{Posterior mean and
#'     variance of the treatment effect (mean difference between treatment and
#'     controls).}
#'   \item{\code{CrI}, \code{Width}}{Credible interval and its length.}
#'   \item{\code{theta_s}}{\code{n_strata}-by-3 matrix of treatment effect
#'     estimates and credible intervals by stratum.}
#'   \item{\code{theta_o}}{1-by-3 matrix of the overall treatment effect and
#'     credible interval.}
#' }
#' When \code{plots = TRUE}, the list also contains the ggplot objects
#' \code{estimate_plot}, \code{ps_group_plot}, \code{ps_group_stratum_plot},
#' \code{density_trace_acf_plot}, \code{overall_mean_plot} and
#' \code{stratum_mean_plot}.
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
#'   \code{\link{bdb_two_arms_bin}}, \code{\link{elastic_fun}},
#'   \code{\link{data_cont}}
#'
#' @examples
#' \donttest{
#' # Requires rjags and JAGS; takes about 15-20 seconds to run.
#' # Trial with N1 = 600 subjects randomized 2:1 (treatment vs control),
#' # borrowing information equivalent to 200 subjects from an external source
#' # of N0 = 3000 patients to achieve a 1:1 ratio of treatment to controls
#' if (requireNamespace("rjags", quietly = TRUE)) {
#'   out <- bdb_two_arms_cont(Y = data_cont$Y,
#'                            Z = data_cont$Z,
#'                            TRT = data_cont$TRT,
#'                            X = data_cont[, paste0('X', 1:10)],
#'                            n_target_borrow = 200,
#'                            n_strata = 3,
#'                            map_a = 5,
#'                            map_b = 1,
#'                            n_ppp = 1000,
#'                            alpha_level = 0.05,
#'                            n_iter = 2000,
#'                            n_adapt = 500,
#'                            n_burn = 500,
#'                            gamma_shape = 0.01,
#'                            gamma_rate = 0.01,
#'                            tau_mean = 0.00001,
#'                            doubleadj = TRUE,
#'                            plots = TRUE,
#'                            seed = 54321)
#'   out$theta_o
#' }
#' }
#' @export
bdb_two_arms_cont <- function(
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
  n_adapt,
  n_burn,
  gamma_shape,
  gamma_rate,
  tau_mean,
  doubleadj,
  plots,
  seed = NULL
) {
  ##############################################################################
  # Function begins
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop(
      "Package 'rjags' is required. Install it and ensure JAGS is installed.",
      call. = FALSE
    )
  }
  # Input checks
  input_checks_cont(
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
    n_adapt,
    n_burn,
    gamma_shape,
    gamma_rate,
    tau_mean,
    doubleadj,
    risk_diff = NULL,
    plots,
    seed
  )

  ##############################################################################
  # Set pseudo-random generator seed
  local_seed(seed)

  ##############################################################################
  # Prior Warning
  if (abs(mean(Y)) > 1e5 || tau_mean * max(mean(Y)^2, safe_var(Y)) > 1) {
    warning(
      "The Normal prior on the mean (centered at 0, precision tau_mean) may be ",
      "informative relative to the scale of Y; consider reducing tau_mean or ",
      "rescaling Y",
      call. = FALSE
    )
  }

  ##############################################################################
  # Define data object
  data <- data.frame(Y, Z, TRT, X)

  ##############################################################################
  # Step 1: Modeling
  ##############################################################################
  # Define exposure mechanism
  mylogit <- fit_ps_model(data, formula = Z ~ . - Y - TRT)

  ##############################################################################
  # Step 2: Propensity Score Calculation
  ##############################################################################
  # Propensity scores
  prd <- ps_predict(mylogit)

  # Add propensity scores to data object
  data <- data.frame(data, prd)

  # Propensity scores of current patients
  prd_current <- prd[Z == 1]

  ##############################################################################
  # Step 3: Trimming
  ##############################################################################
  # Define data_trim and restrict propensity scores
  # to the range of values of subjects in the current trial
  data_trim <- trimming(data, prd, prd_current, sort_by_z = TRUE)

  # Indices current subjects
  is_current <- data_trim$Z == 1
  # Indices external subjects
  is_external <- data_trim$Z == 0
  # Indices current treated subjects
  is_treated <- data_trim$Z == 1 & data_trim$TRT == 1
  # Indices current control subjects
  is_int_control <- data_trim$Z == 1 & data_trim$TRT == 0

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
  # Step 4: Stratification
  ##############################################################################
  # Strata and PS breaks
  strata <- stratify_ps(data_trim$prd, prd1, n_strata)
  data_trim$stratum <- strata$stratum
  q <- strata$breaks

  # n_s0, number of subject in external control by stratum
  n_s0 <- stratum_counts(data_trim$stratum, is_external, n_strata)
  # n_s1, number subjects in the current trial by stratum
  n_s1 <- stratum_counts(data_trim$stratum, is_current, n_strata)
  # n_s11, number subjects in the current trial and treatment arm by stratum
  n_s11 <- stratum_counts(data_trim$stratum, is_treated, n_strata)
  # n_s10, number subjects in the current trial and control arm by stratum
  n_s10 <- stratum_counts(data_trim$stratum, is_int_control, n_strata)

  # Weights for combining results from strata (through a weighted average of trial subjects)
  w_s <- n_s1 / sum(n_s1)

  # There must be at least two subjects on each cell otherwise calculations
  # cannot be computed
  if (any(n_s0 < 2) || any(n_s10 < 2) || any(n_s11 < 2)) {
    stop(
      "There must be at least two subjects on each stratum in the treatment arm, the ",
      "internal control arm and the external data. Consider reducing n_strata or ",
      "revising the propensity score model.",
      call. = FALSE
    )
  }

  # Observed means and 100(1-alpha_level)% confidence intervals
  # Internal Control Arm (ICA)
  r_ica <- mean_CI(data_trim$Y[is_int_control], alpha_level)
  # External Control Arm (ECA)
  r_eca <- mean_CI(data_trim$Y[is_external], alpha_level)
  # ECA by stratum
  r_s0 <- matrix(data = NA, nrow = n_strata, ncol = 3)
  colnames(r_s0) <- c("mean", "lower", "upper")
  for (i in 1:n_strata) {
    r_s0[i, ] <- mean_CI(
      data_trim$Y[is_external & data_trim$stratum == i],
      alpha_level
    )
  }

  # ICA by stratum
  r_s10 <- matrix(data = NA, nrow = n_strata, ncol = 3)
  colnames(r_s10) <- c("mean", "lower", "upper")
  for (i in 1:n_strata) {
    r_s10[i, ] <- mean_CI(
      data_trim$Y[is_int_control & data_trim$stratum == i],
      alpha_level
    )
  }

  ##############################################################################
  # Step 5: Overlapping
  ##############################################################################
  # Calculate vector of overlapping coefficients
  v <- overlap_by_stratum(
    data_trim$prd,
    data_trim$Z,
    data_trim$stratum,
    n_strata
  )

  ##############################################################################
  # Step 6; Weighting
  ##############################################################################
  # Initialize vector of normalized overlapping coefficients
  r <- v / sum(v)

  ##############################################################################
  # Step 7: Calibration by stratum
  ##############################################################################
  # Estimate predictive posterior probability by stratum
  e.m <- rep(0, n_strata) # observed mean in internal control by stratum
  for (i in 1:n_strata) {
    e.m[i] <- mean(data_trim$Y[is_int_control & data_trim$stratum == i])
  }

  PPP <- rep(0, n_strata)
  for (i in 1:n_strata) {
    # Initial values of parameters
    inits_PPP <- list(m = e.m[i])
    # Data list as input to JAGS
    idx_0s <- is_external & data_trim$stratum == i
    idx_10s <- is_int_control & data_trim$stratum == i
    n0_s <- sum(idx_0s)
    y0_s <- data_trim$Y[idx_0s]
    N1_s <- sum(idx_10s)
    data_PPP <- list(
      "n" = n0_s,
      "y" = y0_s,
      "N1" = N1_s,
      "alpha" = gamma_shape,
      "beta" = gamma_rate,
      "tau" = tau_mean
    )
    # Run MCMC by JAGS (predictive mean of N1_s new observations under the external posterior)
    draws_PPP <- run_jags(
      model_string_PPP,
      data = data_PPP,
      inits = inits_PPP,
      params = c("m", "new"),
      n_adapt = n_adapt,
      n_burn = n_burn,
      n_iter = n_ppp,
      seed = if (is.null(seed)) NULL else seed + i
    )
    # One sided measure
    p_one <- mean(draws_PPP$new <= e.m[i])
    PPP[i] <- p_one
  }

  ##############################################################################
  # Step 8: Mapping
  ##############################################################################
  # Mapped predictive posterior probability
  gPPP <- elastic_fun(PPP, map_a, map_b)

  # In case of baseline adjustment only
  if (!doubleadj) {
    gPPP <- rep(1, n_strata)
  }

  ##############################################################################
  # Step 9: Discounting
  ##############################################################################
  # Define baseline discount parameter
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

  # Initialize vector of discount values (current trial: 1, external control: gamma_base*gPPP )
  data_trim$gamma_discount <- 1

  ##############################################################################
  # Steps 10-11 Analysis/Summary
  ##############################################################################
  # Initialize mu_TREAT, mu_CONT, and theta objects
  theta <- matrix(0, nrow = n_strata, ncol = n_iter)
  theta_s <- matrix(0, nrow = n_strata, ncol = 3)
  colnames(theta_s) <- c("mean", "lower", "upper")
  # Obtain posterior mean samples by strata
  for (i in 1:n_strata) {
    # Treatment arm
    # Initial values
    idx_TREAT <- is_treated & data_trim$stratum == i
    mu_ini <- mean(data_trim$Y[idx_TREAT])
    pre_ini <- 1 / safe_var(data_trim$Y[idx_TREAT])
    inits_TREAT <- list(mu = mu_ini, pre = pre_ini)
    # Data list as input to JAGS
    n11_s <- sum(idx_TREAT)
    y11_s <- data_trim$Y[idx_TREAT]
    data_TREAT <- list(
      "n" = n11_s,
      "y" = y11_s,
      "alpha" = gamma_shape,
      "beta" = gamma_rate,
      "tau" = tau_mean
    )
    # Run MCMC by JAGS
    draws_TREAT <- run_jags(
      model_string_TREATMENT,
      data = data_TREAT,
      inits = inits_TREAT,
      params = c("mu", "pre"),
      n_adapt = n_adapt,
      n_burn = n_burn,
      n_iter = n_iter,
      seed = if (is.null(seed)) NULL else seed + 12345 + i
    )
    # Posterior samples
    mu_TREAT_i <- draws_TREAT$mu

    # Internal and external controls
    # Data list
    idx_CONT <- data_trim$TRT == 0 & data_trim$stratum == i
    n_ <- sum(idx_CONT)
    y_ <- data_trim$Y[idx_CONT]
    z_ <- data_trim$Z[idx_CONT]
    gamma_ <- ifelse(z_ == 1, 1, gamma_base_out[i])
    zeros_ <- rep(0, n_)
    data_trim$gamma_discount[idx_CONT] <- gamma_
    data_CONT <- list(
      n = n_,
      y = y_,
      gamma = gamma_,
      zeros = zeros_,
      "alpha" = gamma_shape,
      "beta" = gamma_rate,
      "tau" = tau_mean
    )
    # Initial values
    mu_ <- mean(y_)
    pre_ <- 1 / safe_var(y_)
    inits_CONT <- list(mu = mu_, pre = pre_)
    # Run MCMC by JAGS
    draws_CONT <- run_jags(
      model_string_CONTROL,
      data = data_CONT,
      inits = inits_CONT,
      params = c("mu", "pre"),
      n_adapt = n_adapt,
      n_burn = n_burn,
      n_iter = n_iter,
      seed = if (is.null(seed)) NULL else seed + 23456 + i
    )
    # Posterior samples
    mu_CONT_i <- draws_CONT$mu

    # Difference of posterior distributions by strata
    theta[i, ] <- mu_TREAT_i - mu_CONT_i
    # Posterior mean by stratum and 100(1-alpha_level)% CrI
    theta_s[i, ] <- draws_summary(theta[i, ], alpha_level)
  }

  # Pool results across strata
  pool_results <- pool_strata(theta, w_s, alpha_level)

  # Function output
  output <- list(
    data_trim = data_trim,
    N0_trim = n_external,
    ps_model = mylogit,
    ps_breaks = q,
    n_s0 = n_s0,
    n_s1 = n_s1,
    n_s11 = n_s11,
    n_s10 = n_s10,
    r_s0 = r_s0,
    r_s10 = r_s10,
    r_ica = r_ica,
    r_eca = r_eca,
    A_borrow = A_borrow,
    A_eff = A_eff,
    v = v,
    r = r,
    gamma_base = gamma_base,
    gamma_base_out = gamma_base_out,
    PPP = PPP,
    gPPP = gPPP,
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
  # Plots
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

    # Plot - Combine density, trace, and ACF for trace plot into one plot
    density_plot <- plot_posterior_density(pool_results$theta_m)
    trace_plot <- plot_trace(pool_results$theta_m)
    acf_trace_plot <- acf_plot(pool_results$theta_m, alpha_level)
    density_trace_acf_plot <- ggpubr::ggarrange(
      density_plot,
      trace_plot,
      acf_trace_plot,
      nrow = 3,
      ncol = 1
    )

    # Plot - Sample mean for ICA and ECA
    overall_mean_plot <- plot_group_overall(
      external = r_eca,
      internal = r_ica,
      internal_label = "ICA",
      alpha_level = alpha_level,
      value_label = "Sample mean",
      title = "Overall sample mean for ICA and ECA"
    )

    # Plot - Forest plot of ICA and ECA sample means by stratum
    stratum_mean_plot <- plot_group_by_stratum(
      external = r_s0,
      internal = r_s10,
      internal_label = "ICA",
      n_strata = n_strata,
      alpha_level = alpha_level,
      value_label = "Sample mean",
      title = "Sample mean for ICA and ECA by stratum"
    )

    # Function output
    output <- c(
      output,
      list(
        estimate_plot = estimate_plot,
        ps_group_plot = ps_group_plot,
        ps_group_stratum_plot = ps_group_stratum_plot,
        density_trace_acf_plot = density_trace_acf_plot,
        overall_mean_plot = overall_mean_plot,
        stratum_mean_plot = stratum_mean_plot
      )
    )
  }
  return(output)
}
