valid_bin <- function(two_arms = FALSE) {
  d <- if (two_arms) {
    sim_two_arms(n0 = 20, n1 = 20)
  } else {
    sim_single_arm(n0 = 20, n1 = 20)
  }
  list(
    Y = d$Y,
    Z = d$Z,
    X = d$X,
    TRT = if (two_arms) d$TRT else NULL,
    n_target_borrow = 10,
    n_strata = 2,
    map_a = 1,
    map_b = 1,
    n_ppp = 100,
    alpha_level = 0.05,
    n_iter = 100,
    alpha_shape = 1,
    beta_shape = 1,
    doubleadj = TRUE,
    risk_diff = if (two_arms) NULL else FALSE,
    plots = FALSE,
    seed = NULL
  )
}

valid_cont <- function(two_arms = FALSE) {
  d <- if (two_arms) {
    sim_two_arms(n0 = 20, n1 = 20, binary = FALSE)
  } else {
    sim_single_arm(n0 = 20, n1 = 20, binary = FALSE)
  }
  list(
    Y = d$Y,
    Z = d$Z,
    X = d$X,
    TRT = if (two_arms) d$TRT else NULL,
    n_target_borrow = 10,
    n_strata = 2,
    map_a = 1,
    map_b = 1,
    n_ppp = 100,
    alpha_level = 0.05,
    n_iter = 100,
    n_adapt = 10,
    n_burn = 10,
    gamma_shape = 0.01,
    gamma_rate = 0.01,
    tau_mean = 1e-5,
    doubleadj = TRUE,
    risk_diff = if (two_arms) NULL else FALSE,
    plots = FALSE,
    seed = NULL
  )
}

check_bin <- function(args, ...) {
  do.call(
    input_checks_bin,
    utils::modifyList(args, list(...), keep.null = TRUE)
  )
}
check_cont <- function(args, ...) {
  do.call(
    input_checks_cont,
    utils::modifyList(args, list(...), keep.null = TRUE)
  )
}

test_that("valid inputs pass the checks", {
  expect_no_error(do.call(input_checks_bin, valid_bin()))
  expect_no_error(do.call(input_checks_bin, valid_bin(TRUE)))
  expect_no_error(do.call(input_checks_cont, valid_cont()))
  expect_no_error(do.call(input_checks_cont, valid_cont(TRUE)))
  expect_no_error(check_bin(valid_bin(), seed = 42))
})

test_that("checks on X", {
  a <- valid_bin()
  expect_error(check_bin(a, X = as.vector(a$X)), "matrix or data frame")
  expect_error(check_bin(a, X = unname(a$X)), "non-empty column names")
  X_dup <- a$X
  colnames(X_dup) <- c("X1", "X1")
  expect_error(check_bin(a, X = X_dup), "duplicate column names")
  X_res <- a$X
  colnames(X_res) <- c("X1", "stratum")
  expect_error(check_bin(a, X = X_res), "reserved column names")
  X_chr <- as.data.frame(a$X)
  X_chr$X2 <- as.character(X_chr$X2)
  expect_error(check_bin(a, X = X_chr), "All elements of X must be numeric")
  X_inf <- a$X
  X_inf[1, 1] <- Inf
  expect_error(check_bin(a, X = X_inf), "All elements of X must be finite")
})

test_that("checks on Y, Z and TRT (binary)", {
  a <- valid_bin()
  expect_error(
    check_bin(a, Y = a$Y[-1]),
    "must have equal number of rows",
    fixed = TRUE
  )
  Y_na <- a$Y
  Y_na[1] <- NA
  expect_error(check_bin(a, Y = Y_na), "Missing values are not allowed")
  expect_error(
    check_bin(a, Y = replace(a$Y, 1, 2)),
    "Y must have only 0s or 1s"
  )
  expect_error(
    check_bin(a, Z = replace(a$Z, 1, 2)),
    "Z must have only 0s or 1s"
  )
  expect_error(
    check_bin(a, Z = rep(0, length(a$Z))),
    "No current-trial observations"
  )
  expect_error(
    check_bin(a, Z = rep(1, length(a$Z))),
    "No external observations"
  )

  b <- valid_bin(TRUE)
  expect_error(
    check_bin(b, TRT = replace(b$TRT, 1, 2)),
    "TRT must have only 0s or 1s"
  )
  expect_error(
    check_bin(b, TRT = replace(b$TRT, 1, 1)),
    "External observations must be control treatment only"
  )
  expect_error(
    check_bin(b, TRT = as.numeric(b$Z == 1) * 0),
    "No treated subjects"
  )
  expect_error(check_bin(b, TRT = as.numeric(b$Z == 1)), "No control subjects")
  expect_error(check_bin(b, TRT = b$TRT[-1]), "must have equal number of rows")
})

test_that("checks on Y (continuous)", {
  a <- valid_cont()
  expect_error(
    check_cont(a, Y = as.character(a$Y)),
    "All elements of Y must be numeric"
  )
  expect_error(
    check_cont(a, Y = replace(a$Y, 1, Inf)),
    "All elements of Y must be finite"
  )
  expect_error(
    check_cont(a, Z = replace(a$Z, 1, 0.5)),
    "Z must have only 0s or 1s"
  )
  b <- valid_cont(TRUE)
  expect_error(
    check_cont(b, TRT = replace(b$TRT, 1, 1)),
    "External observations must be control treatment only"
  )
})

test_that("checks on X, Y, Z (continuous, single arm)", {
  a <- valid_cont()
  expect_error(check_cont(a, X = as.vector(a$X)), "matrix or data frame")
  expect_error(check_cont(a, X = unname(a$X)), "non-empty column names")
  expect_error(
    check_cont(a, Y = a$Y[-1]),
    "Y, Z and X must have equal number of rows"
  )
  expect_error(
    check_cont(a, Y = replace(a$Y, 1, NA)),
    "Missing values are not allowed in Y, Z or X"
  )
  X_dup <- a$X
  colnames(X_dup) <- c("X1", "X1")
  expect_error(check_cont(a, X = X_dup), "duplicate column names")
  X_res <- a$X
  colnames(X_res) <- c("X1", "gamma_discount")
  expect_error(check_cont(a, X = X_res), "reserved column names")
  X_chr <- as.data.frame(a$X)
  X_chr$X2 <- as.character(X_chr$X2)
  expect_error(check_cont(a, X = X_chr), "All elements of X must be numeric")
  X_inf <- a$X
  X_inf[1, 1] <- Inf
  expect_error(check_cont(a, X = X_inf), "All elements of X must be finite")
  expect_error(
    check_cont(a, Z = as.character(a$Z)),
    "All elements of Z must be numeric"
  )
  expect_error(
    check_cont(a, Z = rep(0, length(a$Z))),
    "No current-trial observations"
  )
  expect_error(
    check_cont(a, Z = rep(1, length(a$Z))),
    "No external observations"
  )
})

test_that("checks on X, Y, Z and TRT (continuous, two arms)", {
  b <- valid_cont(TRUE)
  expect_error(
    check_cont(b, TRT = b$TRT[-1]),
    "Y, Z, TRT and X must have equal number of rows"
  )
  expect_error(
    check_cont(b, TRT = replace(b$TRT, 1, NA)),
    "Missing values are not allowed in Y, Z, TRT or X"
  )
  X_dup <- b$X
  colnames(X_dup) <- c("X1", "X1")
  expect_error(check_cont(b, X = X_dup), "duplicate column names")
  X_res <- b$X
  colnames(X_res) <- c("X1", "TRT")
  expect_error(check_cont(b, X = X_res), "reserved column names")
  expect_error(
    check_cont(b, Y = as.character(b$Y)),
    "All elements of Y must be numeric"
  )
  expect_error(
    check_cont(b, Y = replace(b$Y, 1, Inf)),
    "All elements of Y must be finite"
  )
  X_chr <- as.data.frame(b$X)
  X_chr$X2 <- as.character(X_chr$X2)
  expect_error(check_cont(b, X = X_chr), "All elements of X must be numeric")
  X_inf <- b$X
  X_inf[1, 1] <- Inf
  expect_error(check_cont(b, X = X_inf), "All elements of X must be finite")
  expect_error(
    check_cont(b, Z = as.character(b$Z)),
    "All elements of Z must be numeric"
  )
  expect_error(
    check_cont(b, TRT = as.character(b$TRT)),
    "All elements of TRT must be numeric"
  )
  expect_error(
    check_cont(b, Z = replace(b$Z, 1, 2)),
    "Z must have only 0s or 1s"
  )
  expect_error(
    check_cont(b, TRT = replace(b$TRT, b$Z == 1, 2)),
    "TRT must have only 0s or 1s"
  )
  expect_error(
    check_cont(b, Z = rep(0, length(b$Z)), TRT = rep(0, length(b$Z))),
    "No current-trial observations"
  )
  expect_error(
    check_cont(b, Z = rep(1, length(b$Z))),
    "No external observations"
  )
  expect_error(check_cont(b, TRT = rep(0, length(b$Z))), "No treated subjects")
  expect_error(check_cont(b, TRT = as.numeric(b$Z == 1)), "No control subjects")
})

test_that("checks on scalar arguments (continuous)", {
  b <- valid_cont()
  expect_error(
    check_cont(b, n_target_borrow = 0),
    "number of subjects to borrow"
  )
  expect_error(
    check_cont(b, n_strata = 0),
    "number of strata must be a positive integer"
  )
  expect_error(
    check_cont(b, map_b = 0),
    "map_b must be a positive finite scalar"
  )
  expect_error(
    check_cont(b, n_ppp = 0.5),
    "supplied n_ppp is not a valid integer"
  )
  expect_error(
    check_cont(b, alpha_level = 1),
    "alpha_level must be in the interval (0,1)",
    fixed = TRUE
  )
  expect_error(
    check_cont(b, n_iter = 0),
    "supplied n_iter is not a valid integer"
  )
  expect_error(check_cont(b, doubleadj = "yes"), "doubleadj must be TRUE/FALSE")
  expect_error(check_cont(b, risk_diff = NA), "risk_diff must be TRUE/FALSE")
  expect_error(check_cont(b, plots = NULL), "plots must be TRUE/FALSE")
  expect_error(
    check_cont(b, seed = -1),
    "seed must be a single non-negative whole number"
  )
})

test_that("checks on scalar tuning arguments", {
  a <- valid_bin()
  msg_borrow <- "number of subjects to borrow must be a positive integer"
  expect_error(check_bin(a, n_target_borrow = 10.5), msg_borrow)
  expect_error(check_bin(a, n_target_borrow = 0), msg_borrow)
  expect_error(check_bin(a, n_target_borrow = c(1, 2)), msg_borrow)
  expect_error(
    check_bin(a, n_strata = 0),
    "number of strata must be a positive integer"
  )
  expect_error(
    check_bin(a, n_strata = 1.5),
    "number of strata must be a positive integer"
  )
  expect_error(
    check_bin(a, map_a = 0),
    "map_a must be a positive finite scalar"
  )
  expect_error(
    check_bin(a, map_b = -1),
    "map_b must be a positive finite scalar"
  )
  expect_error(
    check_bin(a, map_a = Inf),
    "map_a must be a positive finite scalar"
  )
  expect_error(
    check_bin(a, map_b = Inf),
    "map_b must be a positive finite scalar"
  )
  expect_error(
    check_cont(valid_cont(), map_a = Inf),
    "map_a must be a positive finite scalar"
  )
  expect_error(check_bin(a, n_ppp = 0), "supplied n_ppp is not a valid integer")
  expect_error(
    check_bin(a, n_ppp = 10.5),
    "supplied n_ppp is not a valid integer"
  )
  for (al in c(0, 1, 1.5, -0.1)) {
    expect_error(
      check_bin(a, alpha_level = al),
      "alpha_level must be in the interval (0,1)",
      fixed = TRUE
    )
  }
  expect_error(
    check_bin(a, n_iter = 1.5),
    "supplied n_iter is not a valid integer"
  )
  expect_error(
    check_bin(a, alpha_shape = 0),
    "alpha_shape must be a positive finite scalar"
  )
  expect_error(
    check_bin(a, beta_shape = Inf),
    "beta_shape must be a positive finite scalar"
  )
  expect_error(check_bin(a, doubleadj = NA), "doubleadj must be TRUE/FALSE")
  expect_error(check_bin(a, risk_diff = "yes"), "risk_diff must be TRUE/FALSE")
  expect_error(check_bin(a, plots = 1), "plots must be TRUE/FALSE")
  msg_seed <- "seed must be a single non-negative whole number below 2^31"
  for (s in list(c(1, 2), 1.5, -7, NA, Inf, "abc")) {
    expect_error(check_bin(a, seed = s), msg_seed, fixed = TRUE)
  }

  b <- valid_cont()
  expect_error(
    check_cont(b, n_adapt = -1),
    "n_adapt must be a non-negative integer"
  )
  expect_error(
    check_cont(b, n_burn = 2.5),
    "n_burn must be a non-negative integer"
  )
  expect_error(
    check_cont(b, gamma_shape = 0),
    "gamma_shape must be a positive finite scalar"
  )
  expect_error(
    check_cont(b, gamma_rate = -1),
    "gamma_rate must be a positive finite scalar"
  )
  expect_error(
    check_cont(b, tau_mean = 0),
    "tau_mean must be a positive finite scalar"
  )
  expect_no_error(check_cont(b, n_adapt = 0, n_burn = 0))
})

test_that("drivers run the input checks before any computation", {
  d <- sim_single_arm()
  expect_error(
    do.call(
      bdb_single_arm_bin,
      bin_args(d, X = unname(d$X), risk_diff = FALSE)
    ),
    "non-empty column names"
  )
  d2 <- sim_two_arms()
  expect_error(
    do.call(bdb_two_arms_bin, bin_args(d2, map_a = 0)),
    "map_a must be a positive finite scalar"
  )
})
