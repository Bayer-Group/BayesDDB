test_that("binom_ci_exact agrees with stats::binom.test", {
  for (n in c(5, 20, 57)) {
    for (x in unique(c(0, 1, floor(n / 3), n - 1, n))) {
      out <- binom_ci_exact(x, n, 0.05)
      ref <- stats::binom.test(x, n)$conf.int
      expect_equal(unname(out[1, "mean"]), x / n)
      expect_equal(
        unname(out[1, c("lower", "upper")]),
        as.numeric(ref),
        tolerance = 1e-12
      )
    }
  }
  expect_equal(
    colnames(binom_ci_exact(3, 10, 0.1)),
    c("mean", "lower", "upper")
  )
  expect_equal(unname(binom_ci_exact(0, 10, 0.05)[1, "lower"]), 0)
  expect_equal(unname(binom_ci_exact(10, 10, 0.05)[1, "upper"]), 1)
  expect_equal(dim(binom_ci_exact(c(1, 2, 3), c(5, 6, 7), 0.05)), c(3, 3))
  expect_error(binom_ci_exact(0, 0, 0.05), "requires n > 0")
})

test_that("mean_CI returns a t-based confidence interval", {
  set.seed(1)
  x <- rnorm(30, 3, 2)
  out <- mean_CI(x, 0.05)
  expect_equal(colnames(out), c("mean", "lower", "upper"))
  expect_equal(unname(out[1, "mean"]), mean(x))
  expect_equal(
    unname(out[1, c("lower", "upper")]),
    as.numeric(stats::t.test(x, conf.level = 0.95)$conf.int)
  )
  expect_error(mean_CI(1, 0.05), "at least two observations")
})

test_that("draws_summary returns posterior mean and quantiles", {
  set.seed(2)
  d <- rnorm(1000)
  out <- draws_summary(d, 0.1)
  expect_equal(colnames(out), c("mean", "lower", "upper"))
  expect_equal(
    unname(out[1, ]),
    c(mean(d), as.numeric(stats::quantile(d, c(0.05, 0.95))))
  )
})

test_that("pool_strata takes the weighted mean across strata", {
  theta <- rbind(c(1, 2, 3, 4), c(5, 6, 7, 8))
  w_s <- c(1, 3)
  out <- pool_strata(theta, w_s, 0.05)
  expect_equal(out$theta_m, as.vector((1 * theta[1, ] + 3 * theta[2, ]) / 4))
  expect_equal(out$post_mean_theta, mean(out$theta_m))
  expect_equal(out$post_var_theta, stats::var(out$theta_m))
  expect_equal(
    out$CrI,
    as.numeric(stats::quantile(out$theta_m, c(0.025, 0.975)))
  )
  expect_equal(out$Width, out$CrI[2] - out$CrI[1])
  expect_equal(colnames(out$theta_o), c("mean", "lower", "upper"))
  expect_equal(unname(out$theta_o[1, ]), c(out$post_mean_theta, out$CrI))
})

test_that("stratify_ps builds quantile strata of the current-trial scores", {
  set.seed(3)
  prd_current <- runif(60, 0.2, 0.8)
  prd_all <- c(prd_current, runif(100, min(prd_current), max(prd_current)))
  out <- stratify_ps(prd_all, prd_current, 3)
  expect_equal(
    out$breaks,
    as.numeric(stats::quantile(prd_current, c(0, 1 / 3, 2 / 3, 1)))
  )
  expect_setequal(unique(out$stratum), 1:3)
  expect_equal(tabulate(out$stratum[1:60], 3), c(20, 20, 20))

  one <- stratify_ps(prd_all, prd_current, 1)
  expect_equal(one$breaks, range(prd_current))
  expect_true(all(one$stratum == 1L))
})

test_that("stratify_ps errors on degenerate propensity scores", {
  expect_error(
    stratify_ps(c(0.1, 0.5, 0.6), c(0.4, 0.5, 0.6, 0.7), 2),
    "Failed to assign"
  )
  expect_error(stratify_ps(rep(0.5, 5), rep(0.5, 5), 1), "identical")
  expect_error(
    stratify_ps(c(rep(0.5, 10), 0.6), c(rep(0.5, 10), 0.6), 3),
    "not unique"
  )
})

test_that("stratum_counts and stratum_events tabulate by stratum", {
  stratum <- c(1, 1, 2, 2, 3, 3, 3)
  keep <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE)
  y <- c(1, 1, 0, 1, 1, 1, 1)
  expect_equal(stratum_counts(stratum, keep, 3), c(1, 2, 2))
  expect_equal(stratum_counts(stratum, keep, 4), c(1, 2, 2, 0))
  expect_equal(stratum_events(y, stratum, keep, 3), c(1, 1, 2))
  expect_error(stratum_events(c(0, 2), c(1, 1), c(TRUE, TRUE), 1))
})

test_that("trimming keeps subjects within the current-trial score range", {
  data <- data.frame(
    Z = c(0, 0, 1, 1, 0, 0),
    prd = c(0.1, 0.3, 0.3, 0.6, 0.6, 0.9)
  )
  out <- trimming(data, data$prd, data$prd[data$Z == 1])
  expect_equal(out$prd, c(0.3, 0.3, 0.6, 0.6))
  sorted <- trimming(data, data$prd, data$prd[data$Z == 1], sort_by_z = TRUE)
  expect_equal(sorted$Z, c(1, 1, 0, 0))
})

test_that("overlap_by_stratum returns coefficients and rejects zero overlap", {
  set.seed(4)
  prd <- c(runif(40, 0.3, 0.7), runif(40, 0.3, 0.7))
  z <- rep(c(1, 0), each = 40)
  stratum <- ifelse(prd < 0.5, 1, 2)
  v <- overlap_by_stratum(prd, z, stratum, 2)
  expect_length(v, 2)
  expect_true(all(v > 0 & v <= 1))
  expect_error(
    overlap_by_stratum(prd, rep(0, 80), stratum, 2),
    "zero or non-finite"
  )
})

test_that("overlap_by_stratum warns when a single stratum cannot be computed", {
  set.seed(4)
  prd <- runif(80, 0.3, 0.7)
  z <- rep(c(1, 0), each = 40)
  stratum <- ifelse(prd < 0.5, 1, 2)
  z[stratum == 2] <- 0
  expect_warning(
    v <- overlap_by_stratum(prd, z, stratum, 2),
    "could not be computed in stratum 2"
  )
  expect_equal(v[2], 0)
  expect_gt(v[1], 0)
})

test_that("ppp_binary returns probabilities in [0, 1]", {
  set.seed(5)
  out <- ppp_binary(
    ne_s0 = c(10, 20),
    n_s0 = c(30, 40),
    ne_ref = c(5, 30),
    n_ref = c(20, 35),
    alpha_shape = 1,
    beta_shape = 1,
    n_ppp = 500
  )
  expect_length(out, 2)
  expect_true(all(out >= 0 & out <= 1))
  # A reference rate far above the external rate gives a PPP close to 1
  expect_gt(out[2], 0.99)
})

test_that("small utilities behave as documented", {
  expect_true(is_wholenumber(3))
  expect_true(is_wholenumber(3 + 1e-12))
  expect_false(is_wholenumber(3.5))
  expect_false(is_wholenumber(c(1, 2)))
  expect_false(is_wholenumber(NA_real_))
  expect_false(is_wholenumber(Inf))
  expect_false(is_wholenumber("3"))

  expect_equal(safe_var(c(1, 2, 3)), 1)
  expect_equal(safe_var(c(2, 2, 2)), 1e-8)
  expect_equal(safe_var(1), 1e-8)

  expect_equal(stratum_labels(3), c("Stratum 1", "Stratum 2", "Stratum 3"))
})

test_that("fit_ps_model reports failures of the propensity score model", {
  expect_error(
    fit_ps_model(data.frame(Z = c(0, 1, 0, 1)), formula = Z ~ W),
    "Propensity score model failed"
  )
  set.seed(6)
  d <- data.frame(
    Z = rep(c(0, 1), 50),
    Y = rbinom(100, 1, 0.5),
    X1 = rnorm(100)
  )
  fit <- fit_ps_model(d, formula = Z ~ . - Y)
  expect_s3_class(fit, "glm")
  expect_named(stats::coef(fit), c("(Intercept)", "X1"))
  expect_length(ps_predict(fit), 100)

  d_col <- d
  d_col$X2 <- d_col$X1
  expect_error(
    fit_ps_model(d_col, formula = Z ~ . - Y),
    "non-finite coefficients"
  )

  set.seed(7)
  d_sep <- data.frame(Z = rep(c(0, 1), 50), Y = 0)
  d_sep$X1 <- d_sep$Z * 3 + rnorm(100, sd = 0.5)
  expect_error(
    suppressWarnings(fit_ps_model(d_sep, formula = Z ~ . - Y)),
    "did not converge"
  )
})

test_that("ps_predict reports failures and extreme propensity scores", {
  expect_error(ps_predict(list()), "Failed to compute propensity scores")
  near_zero <- stats::lm(y ~ 1, data = data.frame(y = rep(1e-8, 3)))
  expect_warning(ps_predict(near_zero), "numerically near 0 or 1")
})
