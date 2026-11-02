tac_names <- c(
  "data_trim",
  "N0_trim",
  "ps_model",
  "ps_breaks",
  "n_s0",
  "n_s1",
  "n_s11",
  "n_s10",
  "r_s0",
  "r_s10",
  "r_ica",
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

test_that("bdb_two_arms_cont returns the documented elements", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_two_arms(binary = FALSE)
  out <- quiet_jags(do.call(bdb_two_arms_cont, cont_args(d)))
  expect_named(out, tac_names)
  expect_equal(dim(out$theta_s), c(2, 3))
  expect_equal(dim(out$theta_o), c(1, 3))
  expect_equal(out$n_s11 + out$n_s10, out$n_s1)
  expect_length(out$theta_m, 500)
  expect_true(all(out$PPP >= 0 & out$PPP <= 1))
})

test_that("bdb_two_arms_cont discounting follows the documented formulas", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_two_arms(binary = FALSE)
  A <- 30
  out <- quiet_jags(do.call(
    bdb_two_arms_cont,
    cont_args(d, n_target_borrow = A)
  ))
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

test_that("bdb_two_arms_cont recovers a known treatment effect", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_two_arms(n1 = 200, binary = FALSE, delta = 5, sd = 2, seed = 3)
  out <- quiet_jags(do.call(bdb_two_arms_cont, cont_args(d)))
  expect_true(out$theta_o[1, "mean"] > 3 && out$theta_o[1, "mean"] < 7)
  expect_true(out$theta_o[1, "lower"] > 0)
})

test_that("bdb_two_arms_cont is reproducible and supports baseline-only", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_two_arms(binary = FALSE)
  a <- quiet_jags(do.call(bdb_two_arms_cont, cont_args(d, seed = 9)))
  b <- quiet_jags(do.call(bdb_two_arms_cont, cont_args(d, seed = 9)))
  expect_identical(a$theta_m, b$theta_m)
  base <- quiet_jags(do.call(
    bdb_two_arms_cont,
    cont_args(d, doubleadj = FALSE)
  ))
  expect_equal(base$gPPP, rep(1, 2))
  expect_equal(base$gamma_base_out, base$gamma_base)
})

test_that("bdb_two_arms_cont returns plots when requested", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  d <- sim_two_arms(binary = FALSE)
  out <- quiet_jags(do.call(bdb_two_arms_cont, cont_args(d, plots = TRUE)))
  plot_names <- c(
    "estimate_plot",
    "ps_group_plot",
    "ps_group_stratum_plot",
    "density_trace_acf_plot",
    "overall_mean_plot",
    "stratum_mean_plot"
  )
  expect_named(out, c(tac_names, plot_names))
  for (nm in plot_names) {
    expect_true(inherits(out[[nm]], "ggplot"), info = nm)
  }
})
