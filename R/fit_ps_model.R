#' Fit the propensity score model
#'
#' @keywords internal
#' @noRd
fit_ps_model <- function(data, formula) {
  # INPUT
  # data: dataframe containing X, Y , Z and TRT
  # formula: Z ~ . - Y (for single arm designs)
  #          Z ~ . - Y - TRT (for two-arm designs)
  # OUTPUT
  # fitted glm object
  mylogit <- tryCatch(
    eval(bquote(stats::glm(.(formula), family = stats::binomial, data = data))),
    error = function(e) {
      stop("Propensity score model failed: ", e$message, call. = FALSE)
    }
  )

  # Check convergence of PS model
  if (!isTRUE(mylogit$converged)) {
    stop("Propensity score model did not converge", call. = FALSE)
  }
  if (isTRUE(mylogit$boundary)) {
    warning(
      "Propensity score model fit is on the boundary; separation may be present",
      call. = FALSE
    )
  }
  if (!all(is.finite(stats::coef(mylogit)))) {
    stop(
      "Propensity score model produced non-finite coefficients",
      call. = FALSE
    )
  }
  mylogit
}
