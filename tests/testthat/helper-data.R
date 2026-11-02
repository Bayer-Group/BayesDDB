sim_single_arm <- function(
  n0 = 200,
  n1 = 80,
  binary = TRUE,
  p0 = 0.4,
  p1 = 0.35,
  mu0 = 50,
  mu1 = 50,
  sd = 5,
  seed = 1
) {
  set.seed(seed)
  X <- rbind(
    matrix(rnorm(n0 * 2, mean = 0.2), ncol = 2),
    matrix(rnorm(n1 * 2, mean = 0), ncol = 2)
  )
  colnames(X) <- c("X1", "X2")
  Z <- c(rep(0, n0), rep(1, n1))
  Y <- if (binary) {
    c(rbinom(n0, 1, p0), rbinom(n1, 1, p1))
  } else {
    c(rnorm(n0, mu0, sd), rnorm(n1, mu1, sd))
  }
  list(Y = Y, Z = Z, X = X)
}

sim_two_arms <- function(
  n0 = 200,
  n1 = 120,
  binary = TRUE,
  p0 = 0.4,
  p_trt = 0.5,
  mu0 = 50,
  delta = 3,
  sd = 5,
  seed = 1
) {
  d <- sim_single_arm(n0, n1, binary = binary, seed = seed)
  TRT <- c(rep(0, n0), rep(c(0, 1), length.out = n1))
  if (binary) {
    d$Y <- c(
      rbinom(n0, 1, p0),
      rbinom(n1, 1, ifelse(TRT[-(1:n0)] == 1, p_trt, p0))
    )
  } else {
    d$Y <- c(rnorm(n0, mu0, sd), mu0 + delta * TRT[-(1:n0)] + rnorm(n1, 0, sd))
  }
  d$TRT <- TRT
  d
}

bin_args <- function(d, ...) {
  utils::modifyList(
    c(
      d,
      list(
        n_target_borrow = 30,
        n_strata = 3,
        map_a = 1,
        map_b = 1,
        n_ppp = 1000,
        alpha_level = 0.05,
        n_iter = 1000,
        alpha_shape = 1,
        beta_shape = 1,
        doubleadj = TRUE,
        plots = FALSE,
        seed = 123
      )
    ),
    list(...)
  )
}

cont_args <- function(d, ...) {
  utils::modifyList(
    c(
      d,
      list(
        n_target_borrow = 30,
        n_strata = 2,
        map_a = 1,
        map_b = 1,
        n_ppp = 500,
        alpha_level = 0.05,
        n_iter = 500,
        n_adapt = 100,
        n_burn = 100,
        gamma_shape = 0.01,
        gamma_rate = 0.01,
        tau_mean = 1e-5,
        doubleadj = TRUE,
        plots = FALSE,
        seed = 123
      )
    ),
    list(...)
  )
}

collect_warnings <- function(expr) {
  state <- new.env(parent = emptyenv())
  state$msgs <- character(0)
  value <- withCallingHandlers(expr, warning = function(w) {
    state$msgs <- c(state$msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = state$msgs)
}

quiet_jags <- function(expr) {
  res <- NULL
  utils::capture.output({
    res <- expr
  })
  res
}
