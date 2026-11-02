# Regression against reference values computed with the original methodology
# scripts in update/ (see dev/make_reference_fixtures.R), using the bundled data
# sets and the argument values of the original usage examples.

reference <- readRDS(test_path("fixtures", "reference.rds"))

pick_ref <- function(out, ref) out[names(ref)]

# With risk_diff = TRUE and doubleadj = TRUE the original scripts switch
# outcome adjustment off only after A_borrow and A_eff are computed, and return
# the unused gPPP. The package switches it off before discounting, so these
# elements are checked against the documented rule instead.
m1_elements <- c("A_borrow", "A_eff", "gPPP")

expect_baseline_only <- function(out) {
  expect_equal(out$gPPP, rep(1, length(out$gamma_base)))
  expect_identical(out$gamma_base_out, out$gamma_base)
  expect_equal(out$A_borrow, out$n_s0 * out$gamma_base)
  expect_equal(out$A_eff, sum(out$A_borrow))
}

sab_ref_args <- list(
  Y = data_bin_sat$Y,
  Z = data_bin_sat$Z,
  X = data_bin_sat[, paste0("X", 1:10)],
  n_target_borrow = 200,
  n_strata = 3,
  map_a = 1,
  map_b = 3,
  n_ppp = 1e4,
  alpha_level = 0.05,
  n_iter = 1e5,
  alpha_shape = 0.5,
  beta_shape = 0.5,
  doubleadj = TRUE,
  risk_diff = FALSE,
  plots = FALSE,
  seed = 12345
)

sac_ref_args <- list(
  Y = data_cont_sat$Y,
  Z = data_cont_sat$Z,
  X = data_cont_sat[, paste0("X", 1:10)],
  n_target_borrow = 50,
  n_strata = 3,
  map_a = 5,
  map_b = 1,
  n_ppp = 1000,
  alpha_level = 0.05,
  n_iter = 1000,
  n_adapt = 200,
  n_burn = 200,
  gamma_shape = 0.01,
  gamma_rate = 0.01,
  tau_mean = 0.00001,
  doubleadj = TRUE,
  risk_diff = FALSE,
  plots = FALSE,
  seed = 54321
)

test_that("bdb_single_arm_bin reproduces the reference results", {
  ref <- reference$single_arm_bin
  out <- suppressWarnings(do.call(bdb_single_arm_bin, sab_ref_args))
  expect_equal(pick_ref(out, ref), ref)
  expect_equal(round(unname(out$theta_o[1, ]), 4), c(0.3811, 0.3254, 0.4377))
})

test_that("bdb_single_arm_bin risk difference matches the reference without outcome adjustment", {
  ref <- reference$single_arm_bin_rd
  res <- collect_warnings(do.call(
    bdb_single_arm_bin,
    modifyList(sab_ref_args, list(risk_diff = TRUE))
  ))
  out <- res$value
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  keep <- setdiff(names(ref), m1_elements)
  expect_equal(pick_ref(out, ref[keep]), ref[keep])
  expect_baseline_only(out)
  expect_equal(round(unname(out$theta_o[1, ]), 4), c(0.0557, -0.0300, 0.1403))
})

test_that("bdb_two_arms_bin reproduces the reference results", {
  ref <- reference$two_arms_bin
  out <- suppressWarnings(bdb_two_arms_bin(
    Y = data_bin$Y,
    Z = data_bin$Z,
    TRT = data_bin$TRT,
    X = data_bin[, paste0("X", 1:10)],
    n_target_borrow = 100,
    n_strata = 3,
    map_a = 1,
    map_b = 1,
    n_ppp = 1e4,
    alpha_level = 0.05,
    n_iter = 1e5,
    alpha_shape = 0.5,
    beta_shape = 0.5,
    doubleadj = TRUE,
    plots = FALSE,
    seed = 12345
  ))
  expect_equal(pick_ref(out, ref), ref)
  expect_equal(round(unname(out$theta_o[1, ]), 4), c(-0.1042, -0.1926, -0.0143))
})

test_that("bdb_single_arm_cont reproduces the reference results", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  ref <- reference$single_arm_cont
  out <- quiet_jags(suppressWarnings(do.call(
    bdb_single_arm_cont,
    sac_ref_args
  )))
  expect_equal(pick_ref(out, ref), ref)
})

test_that("bdb_single_arm_cont mean difference matches the reference without outcome adjustment", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  ref <- reference$single_arm_cont_rd
  res <- collect_warnings(quiet_jags(do.call(
    bdb_single_arm_cont,
    modifyList(sac_ref_args, list(risk_diff = TRUE))
  )))
  out <- res$value
  expect_true(any(grepl("Double adjustment is not valid", res$warnings)))
  keep <- setdiff(names(ref), m1_elements)
  expect_equal(pick_ref(out, ref[keep]), ref[keep])
  expect_baseline_only(out)
})

test_that("bdb_two_arms_cont reproduces the reference results", {
  skip_if_not_installed("rjags")
  skip_on_cran()
  ref <- reference$two_arms_cont
  out <- quiet_jags(suppressWarnings(bdb_two_arms_cont(
    Y = data_cont$Y,
    Z = data_cont$Z,
    TRT = data_cont$TRT,
    X = data_cont[, paste0("X", 1:10)],
    n_target_borrow = 200,
    n_strata = 3,
    map_a = 5,
    map_b = 1,
    n_ppp = 1000,
    alpha_level = 0.05,
    n_iter = 1000,
    n_adapt = 200,
    n_burn = 200,
    gamma_shape = 0.01,
    gamma_rate = 0.01,
    tau_mean = 0.00001,
    doubleadj = TRUE,
    plots = FALSE,
    seed = 54321
  )))
  expect_equal(pick_ref(out, ref), ref)
})
