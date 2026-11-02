sab_names <- c(
  "data_trim",
  "N0_trim",
  "ps_model",
  "ps_breaks",
  "n_s1",
  "n_s0",
  "ne_s1",
  "ne_s0",
  "r_s1",
  "r_s0",
  "r_sat",
  "r_eca",
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

test_that("bdb_single_arm_bin returns the documented elements", {
  d <- sim_single_arm()
  out <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE))
  expect_type(out, "list")
  expect_named(out, sab_names)
  expect_equal(dim(out$theta_s), c(3, 3))
  expect_equal(dim(out$theta_o), c(1, 3))
  expect_equal(colnames(out$theta_o), c("mean", "lower", "upper"))
  expect_equal(dim(out$r_s1), c(3, 3))
  expect_length(out$theta_m, 1000)
  expect_length(out$ps_breaks, 4)
  expect_s3_class(out$ps_model, "glm")
  expect_true(all(
    c("prd", "stratum", "gamma_discount") %in% names(out$data_trim)
  ))
  expect_true(out$theta_o[1, "lower"] <= out$theta_o[1, "mean"])
  expect_true(out$theta_o[1, "upper"] >= out$theta_o[1, "mean"])
  expect_true(all(out$PPP >= 0 & out$PPP <= 1))
  expect_equal(out$N0_trim, sum(out$n_s0))
})

test_that("bdb_single_arm_bin discounting follows the documented formulas", {
  d <- sim_single_arm()
  A <- 30
  out <- do.call(
    bdb_single_arm_bin,
    bin_args(d, n_target_borrow = A, risk_diff = FALSE)
  )
  expect_equal(out$r, out$v / sum(out$v))
  expect_equal(out$gPPP, elastic_fun(out$PPP, 1, 1))
  expect_equal(out$gamma_base, pmin(1, A * out$r / out$n_s0))
  expect_equal(out$gamma_base_out, pmin(1, A * out$r * out$gPPP / out$n_s0))
  expect_equal(out$A_borrow, out$n_s0 * out$gamma_base_out)
  expect_equal(out$A_eff, sum(out$A_borrow))
  expect_equal(out$w_s, out$n_s1 / sum(out$n_s1))
  expect_equal(unname(out$theta_o[1, "mean"]), out$post_mean_theta)
  expect_equal(out$Width, out$CrI[2] - out$CrI[1])

  dt <- out$data_trim
  expect_true(all(dt$gamma_discount[dt$Z == 1] == 1))
  ext <- dt$Z == 0
  expect_equal(dt$gamma_discount[ext], out$gamma_base_out[dt$stratum[ext]])
})

test_that("bdb_single_arm_bin with doubleadj = FALSE adjusts for baseline only", {
  d <- sim_single_arm()
  out <- do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = FALSE, risk_diff = FALSE)
  )
  expect_equal(out$gPPP, rep(1, 3))
  expect_equal(out$gamma_base_out, out$gamma_base)
})

test_that("bdb_single_arm_bin is reproducible with a seed", {
  d <- sim_single_arm()
  a <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 11))
  b <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 11))
  c <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 12))
  expect_identical(a$theta_m, b$theta_m)
  expect_identical(a$PPP, b$PPP)
  expect_false(identical(a$theta_m, c$theta_m))
  expect_no_error(do.call(
    bdb_single_arm_bin,
    bin_args(d, risk_diff = FALSE, seed = NULL)
  ))
})

test_that("bdb_single_arm_bin risk_diff changes only the estimand", {
  d <- sim_single_arm()
  pooled <- do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = FALSE, risk_diff = FALSE)
  )
  rdiff <- do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = FALSE, risk_diff = TRUE)
  )
  expect_identical(pooled$PPP, rdiff$PPP)
  expect_identical(pooled$gamma_base_out, rdiff$gamma_base_out)
  expect_false(isTRUE(all.equal(pooled$theta_o, rdiff$theta_o)))
  # A pooled rate lies in [0, 1]; the risk difference between similar arms is near 0
  expect_true(pooled$theta_o[1, "mean"] > 0 && pooled$theta_o[1, "mean"] < 1)
  expect_true(abs(rdiff$theta_o[1, "mean"]) < 0.3)
})

test_that("bdb_single_arm_bin risk_diff does not apply outcome adjustment", {
  d <- sim_single_arm()
  res <- collect_warnings(do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = TRUE, risk_diff = TRUE)
  ))
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  out <- res$value
  expect_equal(out$gPPP, rep(1, 3))
  expect_identical(out$gamma_base_out, out$gamma_base)
  expect_equal(out$A_borrow, out$n_s0 * out$gamma_base)
  expect_equal(out$A_eff, sum(out$A_borrow))
  dt <- out$data_trim
  ext <- dt$Z == 0
  expect_equal(dt$gamma_discount[ext], out$gamma_base[dt$stratum[ext]])

  # Identical to requesting baseline adjustment only
  base <- do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = FALSE, risk_diff = TRUE)
  )
  expect_identical(out$PPP, base$PPP)
  expect_identical(out$theta_m, base$theta_m)
  expect_identical(out$theta_o, base$theta_o)

  # The pooled estimand still uses outcome adjustment and does not warn
  pooled <- collect_warnings(do.call(
    bdb_single_arm_bin,
    bin_args(d, doubleadj = TRUE, risk_diff = FALSE)
  ))
  expect_false(any(grepl("Double adjustment is not valid", pooled$warnings)))
  expect_equal(pooled$value$gPPP, elastic_fun(pooled$value$PPP, 1, 1))
})

test_that("bdb_single_arm_bin risk difference is not attenuated by a treatment effect", {
  # True risk difference 0.25: the outcome difference is the treatment effect
  # and must not switch borrowing off
  d <- sim_single_arm(n0 = 400, n1 = 120, p0 = 0.35, p1 = 0.6)
  out <- suppressWarnings(do.call(
    bdb_single_arm_bin,
    bin_args(d, n_target_borrow = 120, risk_diff = TRUE)
  ))
  expect_equal(out$A_eff, sum(out$n_s0 * out$gamma_base))
  expect_true(out$A_eff > 60)
  expect_true(out$theta_o[1, "lower"] < 0.25 && out$theta_o[1, "upper"] > 0.25)
  expect_true(out$theta_o[1, "lower"] > 0)
})

test_that("bdb_single_arm_bin works with a single stratum", {
  d <- sim_single_arm()
  out <- do.call(
    bdb_single_arm_bin,
    bin_args(d, n_strata = 1, risk_diff = FALSE)
  )
  expect_equal(dim(out$theta_s), c(1, 3))
  expect_equal(out$w_s, 1)
  expect_true(all(out$data_trim$stratum == 1))
})

test_that("bdb_single_arm_bin warns and errors as documented", {
  d <- sim_single_arm()
  res <- collect_warnings(
    do.call(
      bdb_single_arm_bin,
      bin_args(d, n_target_borrow = 10000, risk_diff = FALSE)
    )
  )
  expect_true(any(grepl(
    "exceeds number of trimmed external subjects",
    res$warnings
  )))
  expect_true(all(res$value$gamma_base == 1))

  d_diff <- sim_single_arm(p0 = 0.9, p1 = 0.1)
  res <- collect_warnings(
    do.call(bdb_single_arm_bin, bin_args(d_diff, risk_diff = TRUE))
  )
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  expect_false(any(grepl("Large outcome difference", res$warnings)))
  expect_equal(res$value$gamma_base_out, res$value$gamma_base)

  expect_error(
    do.call(
      bdb_single_arm_bin,
      bin_args(sim_single_arm(n1 = 30), n_strata = 20, risk_diff = FALSE)
    ),
    "at least two subjects on each stratum"
  )
})

test_that("bdb_single_arm_bin returns plots when requested", {
  d <- sim_single_arm()
  out <- do.call(
    bdb_single_arm_bin,
    bin_args(d, risk_diff = FALSE, plots = TRUE)
  )
  plot_names <- c(
    "estimate_plot",
    "ps_group_plot",
    "ps_group_stratum_plot",
    "density_plot",
    "overall_mean_plot",
    "stratum_mean_plot"
  )
  expect_named(out, c(sab_names, plot_names))
  for (nm in plot_names) {
    expect_true(inherits(out[[nm]], "ggplot"), info = nm)
  }
})

test_that("the seed argument makes results reproducible and restores the caller's RNG state", {
  d <- sim_single_arm()
  set.seed(2024)
  expected_next <- stats::runif(1)
  set.seed(2024)
  out1 <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 7))
  expect_equal(stats::runif(1), expected_next)
  out2 <- do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 7))
  expect_identical(out1$theta_m, out2$theta_m)

  old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  rm(".Random.seed", envir = globalenv())
  do.call(bdb_single_arm_bin, bin_args(d, risk_diff = FALSE, seed = 7))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})
