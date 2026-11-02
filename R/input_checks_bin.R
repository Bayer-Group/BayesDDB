#' Input checks for the binary-outcome functions
#'
#' @keywords internal
#' @noRd
input_checks_bin <- function(
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
  risk_diff = NULL,
  plots,
  seed
) {
  # Check for matrix and data frame
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data frame of confounders")
  }

  # Check for null or empty names in X
  if (is.null(colnames(X)) || any(colnames(X) == "")) {
    stop("X must have non-empty column names")
  }

  # Convert matrix X into a data frame
  X <- as.data.frame(X)

  # Single-arm branch
  if (is.null(TRT)) {
    # Check X, Y, Z
    if (length(unique(c(length(Y), length(Z), nrow(X)))) != 1) {
      stop("Y, Z and X must have equal number of rows")
    }

    if (anyNA(Y) || anyNA(Z) || anyNA(X)) {
      stop("Missing values are not allowed in Y, Z or X")
    }

    if (ncol(X) < 1) {
      stop("X must have at least one confounder")
    }

    if (anyDuplicated(colnames(X))) {
      stop("X contains duplicate column names")
    }

    if (
      any(
        colnames(X) %in% c("Y", "Z", "TRT", "prd", "stratum", "gamma_discount")
      )
    ) {
      stop("X contains reserved column names")
    }

    if (!is.numeric(Y)) {
      stop("All elements of Y must be numeric")
    }

    if (!all(is.finite(Y))) {
      stop("All elements of Y must be finite")
    }

    if (!all(sapply(X, is.numeric))) {
      stop("All elements of X must be numeric")
    }

    if (!all(sapply(X, is.finite))) {
      stop("All elements of X must be finite")
    }

    if (!is.numeric(Z)) {
      stop("All elements of Z must be numeric")
    }

    if (!all(Y == 0 | Y == 1)) {
      stop("Y must have only 0s or 1s")
    }

    if (!all(Z == 0 | Z == 1)) {
      stop("Z must have only 0s or 1s")
    }

    if (sum(Z == 1) == 0) {
      stop("No current-trial observations found")
    }

    if (sum(Z == 0) == 0) stop("No external observations found")
  } else {
    # Two-arms branch

    # Check X, Y, Z
    if (length(unique(c(length(Y), length(Z), length(TRT), nrow(X)))) != 1) {
      stop("Y, Z, TRT and X must have equal number of rows")
    }

    if (anyNA(Y) || anyNA(Z) || anyNA(TRT) || anyNA(X)) {
      stop("Missing values are not allowed in Y, Z, TRT or X")
    }

    if (ncol(X) < 1) {
      stop("X must have at least one confounder")
    }

    if (is.null(colnames(X)) || any(colnames(X) == "")) {
      stop("X must have non-empty column names")
    }

    if (anyDuplicated(colnames(X))) {
      stop("X contains duplicate column names")
    }

    if (
      any(
        colnames(X) %in% c("Y", "Z", "TRT", "prd", "stratum", "gamma_discount")
      )
    ) {
      stop("X contains reserved column names")
    }

    if (!is.numeric(Y)) {
      stop("All elements of Y must be numeric")
    }

    if (!all(is.finite(Y))) {
      stop("All elements of Y must be finite")
    }

    if (!all(sapply(X, is.numeric))) {
      stop("All elements of X must be numeric")
    }

    if (!all(sapply(X, is.finite))) {
      stop("All elements of X must be finite")
    }

    if (!is.numeric(Z)) {
      stop("All elements of Z must be numeric")
    }

    if (!is.numeric(TRT)) {
      stop("All elements of TRT must be numeric")
    }

    if (!all(Y == 0 | Y == 1)) {
      stop("Y must have only 0s or 1s")
    }

    if (!all(Z == 0 | Z == 1)) {
      stop("Z must have only 0s or 1s")
    }

    if (!all(TRT == 0 | TRT == 1)) {
      stop("TRT must have only 0s or 1s")
    }

    if (any(Z == 0 & TRT != 0)) {
      stop("External observations must be control treatment only")
    }

    if (sum(Z == 1) == 0) {
      stop("No current-trial observations found")
    }

    if (sum(Z == 0) == 0) {
      stop("No external observations found")
    }

    if (sum(Z == 1 & TRT == 1) == 0) {
      stop("No treated subjects in current trial")
    }

    if (sum(Z == 1 & TRT == 0) == 0) {
      stop("No control subjects in current trial")
    }
  }

  # Check n_target_borrow
  if (
    length(n_target_borrow) != 1 ||
      !is.numeric(n_target_borrow) ||
      is.na(n_target_borrow) ||
      n_target_borrow <= 0 ||
      !is_wholenumber(n_target_borrow)
  ) {
    stop("The number of subjects to borrow must be a positive integer")
  }

  # Check n_strata
  if (
    length(n_strata) != 1 ||
      !is.numeric(n_strata) ||
      is.na(n_strata) ||
      n_strata < 1 ||
      !is_wholenumber(n_strata)
  ) {
    stop("The number of strata must be a positive integer")
  }

  # Check tuning parameters map_a and map_b
  if (
    !is.numeric(map_a) ||
      length(map_a) != 1 ||
      is.na(map_a) ||
      !is.finite(map_a) ||
      map_a <= 0
  ) {
    stop("map_a must be a positive finite scalar")
  }

  if (
    !is.numeric(map_b) ||
      length(map_b) != 1 ||
      is.na(map_b) ||
      !is.finite(map_b) ||
      map_b <= 0
  ) {
    stop("map_b must be a positive finite scalar")
  }

  # Check n_ppp
  if (
    !is.numeric(n_ppp) ||
      !is_wholenumber(n_ppp) ||
      length(n_ppp) != 1 ||
      is.na(n_ppp) ||
      !is.finite(n_ppp) ||
      n_ppp < 1 ||
      n_ppp > 2^31 - 1
  ) {
    stop("supplied n_ppp is not a valid integer")
  }

  # Check alpha_level
  if (
    !is.numeric(alpha_level) ||
      length(alpha_level) != 1 ||
      is.na(alpha_level) ||
      alpha_level <= 0 ||
      alpha_level >= 1
  ) {
    stop("alpha_level must be in the interval (0,1)")
  }

  # Check n_iter
  if (
    !is.numeric(n_iter) ||
      !is_wholenumber(n_iter) ||
      length(n_iter) != 1 ||
      is.na(n_iter) ||
      !is.finite(n_iter) ||
      n_iter < 1 ||
      n_iter > 2^31 - 1
  ) {
    stop("supplied n_iter is not a valid integer")
  }

  # Check alpha shape parameters for prior
  if (
    !is.numeric(alpha_shape) ||
      length(alpha_shape) != 1 ||
      is.na(alpha_shape) ||
      !is.finite(alpha_shape) ||
      alpha_shape <= 0
  ) {
    stop("alpha_shape must be a positive finite scalar")
  }

  # Check beta shape parameters for prior
  if (
    !is.numeric(beta_shape) ||
      length(beta_shape) != 1 ||
      is.na(beta_shape) ||
      !is.finite(beta_shape) ||
      beta_shape <= 0
  ) {
    stop("beta_shape must be a positive finite scalar")
  }

  # Check doubleadj
  if (!is.logical(doubleadj) || length(doubleadj) != 1 || is.na(doubleadj)) {
    stop("doubleadj must be TRUE/FALSE")
  }

  # Check risk_diff
  if (!is.null(risk_diff)) {
    if (!is.logical(risk_diff) || length(risk_diff) != 1 || is.na(risk_diff)) {
      stop("risk_diff must be TRUE/FALSE")
    }
  }

  # Check plots
  if (!is.logical(plots) || length(plots) != 1 || is.na(plots)) {
    stop("plots must be TRUE/FALSE")
  }

  # Check seed
  if (!is.null(seed)) {
    if (
      !is.numeric(seed) ||
        !is_wholenumber(seed) ||
        length(seed) != 1 ||
        is.na(seed) ||
        !is.finite(seed) ||
        seed < 0 ||
        seed > 2^31 - 1
    ) {
      stop("seed must be a single non-negative whole number below 2^31")
    }
  }
}
