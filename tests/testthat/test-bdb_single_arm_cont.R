sac_names <- c(
  "data_trim",
  "N0_trim",
  "ps_model",
  "ps_breaks",
  "n_s0",
  "n_s1",
  "r_s0",
  "r_s1",
  "r_sat",
  "r_eca",
  "A_borrow",
  "A_eff",
  "v",
  "r",
  "gamma_base",
  "gamma_base_out",
  "PPP",
  "gPPP",
  "w_s",
  "theta_m",
  "post_mean_theta",
  "post_var_theta",
  "CrI",
  "Width",
  "theta_s",
  "theta_o"
)

test_that("bdb_single_arm_cont returns the documented elements", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  out <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE)
  ))
  expect_named(out, sac_names)
  expect_equal(dim(out$theta_s), c(2, 3))
  expect_equal(dim(out$theta_o), c(1, 3))
  expect_length(out$theta_m, 500)
  expect_true(all(out$PPP >= 0 & out$PPP <= 1))
  # The pooled mean sits on the data scale
  expect_true(out$theta_o[1, "mean"] > 45 && out$theta_o[1, "mean"] < 55)
})

test_that("bdb_single_arm_cont discounting follows the documented formulas", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  A <- 30
  out <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, n_target_borrow = A, risk_diff = FALSE)
  ))
  expect_equal(out$r, out$v / sum(out$v))
  expect_equal(out$gPPP, elastic_fun(out$PPP, 1, 1))
  expect_equal(out$gamma_base, pmin(1, A * out$r / out$n_s0))
  expect_equal(out$gamma_base_out, pmin(1, A * out$r * out$gPPP / out$n_s0))
  expect_equal(out$A_eff, sum(out$n_s0 * out$gamma_base_out))
  expect_equal(out$w_s, out$n_s1 / sum(out$n_s1))
  dt <- out$data_trim
  ext <- dt$Z == 0
  expect_true(all(dt$gamma_discount[!ext] == 1))
  expect_equal(dt$gamma_discount[ext], out$gamma_base_out[dt$stratum[ext]])
})

test_that("bdb_single_arm_cont risk_diff returns a mean difference", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  pooled <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, doubleadj = FALSE, risk_diff = FALSE)
  ))
  rdiff <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, doubleadj = FALSE, risk_diff = TRUE)
  ))
  expect_identical(pooled$PPP, rdiff$PPP)
  expect_identical(pooled$gamma_base_out, rdiff$gamma_base_out)
  expect_true(abs(rdiff$theta_o[1, "mean"]) < 5)
  expect_true(rdiff$theta_o[1, "lower"] <= rdiff$theta_o[1, "mean"])
  expect_true(rdiff$theta_o[1, "upper"] >= rdiff$theta_o[1, "mean"])
})

test_that("bdb_single_arm_cont risk_diff does not apply outcome adjustment", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  res <- collect_warnings(quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, doubleadj = TRUE, risk_diff = TRUE)
  )))
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  out <- res$value
  expect_equal(out$gPPP, rep(1, 2))
  expect_identical(out$gamma_base_out, out$gamma_base)
  expect_equal(out$A_eff, sum(out$n_s0 * out$gamma_base))

  base <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, doubleadj = FALSE, risk_diff = TRUE)
  ))
  expect_identical(out$PPP, base$PPP)
  expect_identical(out$theta_m, base$theta_m)
})

test_that("bdb_single_arm_cont is reproducible and supports baseline-only", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  a <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE, seed = 3)
  ))
  b <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE, seed = 3)
  ))
  expect_identical(a$theta_m, b$theta_m)
  expect_identical(a$PPP, b$PPP)

  base <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE, doubleadj = FALSE)
  ))
  expect_equal(base$gPPP, rep(1, 2))
  expect_equal(base$gamma_base_out, base$gamma_base)
})

test_that("bdb_single_arm_cont returns plots when requested", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  out <- quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE, plots = TRUE)
  ))
  plot_names <- c(
    "estimate_plot",
    "ps_group_plot",
    "ps_group_stratum_plot",
    "density_trace_acf_plot",
    "overall_mean_plot",
    "stratum_mean_plot"
  )
  expect_named(out, c(sac_names, plot_names))
  for (nm in plot_names) {
    expect_true(inherits(out[[nm]], "ggplot"), info = nm)
  }
})

test_that("continuous functions accept seeds up to the documented upper bound", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  big_seed <- 2^31 - 1
  d <- sim_single_arm(binary = FALSE)
  expect_no_warning(quiet_jags(do.call(
    bdb_single_arm_cont,
    cont_args(d, risk_diff = FALSE, seed = big_seed)
  )))
  d2 <- sim_two_arms(binary = FALSE)
  expect_no_warning(quiet_jags(do.call(
    bdb_two_arms_cont,
    cont_args(d2, seed = big_seed)
  )))
})

test_that("continuous functions print nothing to the console", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  expect_silent(do.call(bdb_single_arm_cont, cont_args(d, risk_diff = FALSE)))
  d2 <- sim_two_arms(binary = FALSE)
  expect_silent(do.call(bdb_two_arms_cont, cont_args(d2)))
})

test_that("bdb_single_arm_cont warns about an informative prior and keeps borrowing under a mean difference", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_single_arm(binary = FALSE)
  expect_warning(
    do.call(bdb_single_arm_cont, cont_args(d, risk_diff = FALSE, tau_mean = 1)),
    "may be informative relative to the scale of Y"
  )
  # A true mean difference of 10 is the treatment effect: borrowing is kept
  d_diff <- sim_single_arm(binary = FALSE, mu0 = 50, mu1 = 60)
  res <- collect_warnings(do.call(
    bdb_single_arm_cont,
    cont_args(d_diff, risk_diff = TRUE)
  ))
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  expect_false(any(grepl("Large outcome difference", res$warnings)))
  expect_true(all(res$value$PPP > 0.9))
  expect_equal(res$value$gamma_base_out, res$value$gamma_base)
  ci <- res$value$theta_o[1, ]
  expect_true(ci["lower"] < 10 && ci["upper"] > 10)
})
