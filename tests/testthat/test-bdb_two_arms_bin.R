tab_names <- c(
  "data_trim",
  "N0_trim",
  "ps_model",
  "ps_breaks",
  "n_s1",
  "n_s0",
  "n_s11",
  "n_s10",
  "ne_s11",
  "ne_s10",
  "ne_s0",
  "r_s11",
  "r_s10",
  "r_s0",
  "r_ica",
  "r_eca",
  "theta_n",
  "v",
  "r",
  "A_borrow",
  "A_eff",
  "PPP",
  "gPPP",
  "gamma_base",
  "gamma_base_out",
  "w_s",
  "theta_m",
  "post_mean_theta",
  "post_var_theta",
  "CrI",
  "Width",
  "theta_s",
  "theta_o"
)

test_that("bdb_two_arms_bin returns the documented elements", {
  d <- sim_two_arms()
  out <- do.call(bdb_two_arms_bin, bin_args(d))
  expect_type(out, "list")
  expect_named(out, tab_names)
  expect_equal(dim(out$theta_s), c(3, 3))
  expect_equal(dim(out$theta_o), c(1, 3))
  expect_named(out$theta_n, c("mean", "lower", "upper"))
  expect_equal(out$n_s11 + out$n_s10, out$n_s1)
  expect_length(out$theta_m, 1000)
  expect_true(all(out$PPP >= 0 & out$PPP <= 1))
  dt <- out$data_trim
  expect_equal(
    out$ne_s11,
    as.vector(tabulate(dt$stratum[dt$Z == 1 & dt$TRT == 1 & dt$Y == 1], 3))
  )
})

test_that("bdb_two_arms_bin discounting follows the documented formulas", {
  d <- sim_two_arms()
  A <- 30
  out <- do.call(
    bdb_two_arms_bin,
    bin_args(d, n_target_borrow = A, map_a = 5, map_b = 3)
  )
  expect_equal(out$r, out$v / sum(out$v))
  expect_equal(out$gPPP, elastic_fun(out$PPP, 5, 3))
  expect_equal(out$gamma_base, pmin(1, A * out$r / out$n_s0))
  expect_equal(out$gamma_base_out, pmin(1, A * out$r * out$gPPP / out$n_s0))
  expect_equal(out$A_eff, sum(out$n_s0 * out$gamma_base_out))
  expect_equal(out$w_s, out$n_s1 / sum(out$n_s1))

  dt <- out$data_trim
  expect_true(all(dt$gamma_discount[dt$Z == 1] == 1))
  ext <- dt$Z == 0
  expect_equal(dt$gamma_discount[ext], out$gamma_base_out[dt$stratum[ext]])
})

test_that("bdb_two_arms_bin with doubleadj = FALSE adjusts for baseline only", {
  d <- sim_two_arms()
  out <- do.call(bdb_two_arms_bin, bin_args(d, doubleadj = FALSE))
  expect_equal(out$gPPP, rep(1, 3))
  expect_equal(out$gamma_base_out, out$gamma_base)
})

test_that("bdb_two_arms_bin is reproducible with a seed", {
  d <- sim_two_arms()
  a <- do.call(bdb_two_arms_bin, bin_args(d, seed = 5))
  b <- do.call(bdb_two_arms_bin, bin_args(d, seed = 5))
  expect_identical(a$theta_m, b$theta_m)
  expect_identical(a$theta_n, b$theta_n)
})

test_that("bdb_two_arms_bin detects a treatment effect", {
  d <- sim_two_arms(n1 = 300, p0 = 0.3, p_trt = 0.7, seed = 2)
  out <- do.call(bdb_two_arms_bin, bin_args(d))
  expect_true(out$theta_o[1, "lower"] > 0)
  expect_true(out$theta_n["lower"] > 0)
})

test_that("bdb_two_arms_bin errors when strata are too small", {
  d <- sim_two_arms(n1 = 30)
  expect_error(
    do.call(bdb_two_arms_bin, bin_args(d, n_strata = 10)),
    "at least two subjects on each stratum"
  )
})

test_that("bdb_two_arms_bin returns plots when requested", {
  d <- sim_two_arms()
  out <- do.call(bdb_two_arms_bin, bin_args(d, plots = TRUE))
  plot_names <- c(
    "estimate_plot",
    "ps_group_plot",
    "ps_group_stratum_plot",
    "density_plot",
    "overall_mean_plot",
    "stratum_mean_plot"
  )
  expect_named(out, c(tab_names, plot_names))
  for (nm in plot_names) {
    expect_true(inherits(out[[nm]], "ggplot"), info = nm)
  }
})
