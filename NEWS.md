# BayesDDB NEWS

## BayesDDB 0.5.0

This release incorporates the source code release of 25.09.2026, 
which addresses the review finding that outcome adjustment is not valid 
for the difference estimand of single-arm designs.

### Breaking changes

- **No outcome adjustment for `risk_diff = TRUE`.** In
  `bdb_single_arm_bin()` and `bdb_single_arm_cont()`, the external data are
  the comparator when `risk_diff = TRUE`, so the posterior predictive
  probability compares the treated trial arm with the external controls and
  measures the treatment effect rather than prior-data conflict. Previously,
  a true treatment effect reduced borrowing and attenuated the estimated
  difference. With `risk_diff = TRUE` and `doubleadj = TRUE`, both functions
  now issue a warning and apply baseline adjustment only: `gPPP` is set to 1,
  `gamma_base_out` equals `gamma_base`, and `A_borrow` and `A_eff` report the
  borrowing actually used. `PPP` is still computed and returned. Results with
  `risk_diff = TRUE` and `doubleadj = TRUE` therefore change; they are now
  identical to those with `doubleadj = FALSE`.
- The warning "Large outcome difference in at least one stratum strongly
  limits borrowing" under `risk_diff = TRUE` has been removed, because the
  situation it described can no longer occur.
- **`data_bin_sat` has been replaced** by the data set of the source code
  release: 2200 rows (200 current-trial and 2000 external subjects) instead
  of 1200 rows. Results of examples based on `data_bin_sat` change.

### Other changes

- In the source code release, outcome adjustment is switched off only after
  `A_borrow` and `A_eff` have been computed, so these elements and `gPPP`
  did not match the discount actually used. The package applies the change
  before the discounting step, so all returned discounting quantities are
  consistent.
- The examples of `bdb_single_arm_bin()` and `bdb_single_arm_cont()` now use
  `risk_diff = FALSE`, as in the usage examples of the source code release.
- Documentation of `doubleadj`, `risk_diff`, `PPP` and `gPPP` and the
  vignette section "Single-Arm Trials: Choosing the Estimand" explain why
  outcome adjustment is not applied for the difference estimand.
- The reference fixtures (`tests/testthat/fixtures/reference.rds`) are
  regenerated from the source code release with `dev/make_reference_fixtures.R`,
  and `dev/equivalence_check.R` now compares the package with that release.
- `RoxygenNote` has been removed from `DESCRIPTION`
  (`Config/roxygen2/version` is kept).

## BayesDDB 0.4.0

This release replaces the internals of the package with a functionally updated
and fundamentally refactored implementation of the methodology. Numerical
results change, and the interface is not backward compatible.

### Breaking changes

- **Arguments renamed** in all analysis functions, following the refactored
  implementation: `A` -> `n_target_borrow`, `s` -> `n_strata`, `a` -> `map_a`,
  `p` -> `map_b`, `za` -> `alpha_level`, `K` -> `n_iter` (binary functions),
  `alpha`/`beta` -> `alpha_shape`/`beta_shape` (binary functions) and
  `gamma_shape`/`gamma_rate` (continuous functions), `tau` -> `tau_mean`,
  `double_adj` -> `doubleadj`, `show_plots` -> `plots`. New arguments are
  `n_ppp`, `seed` and, for the continuous functions, `n_adapt` and `n_burn`.
- **No default values** apart from `n_ppp = 1e3` and `seed = NULL`; all other
  arguments must be supplied.
- `elastic_function(ppp, a, p)` is replaced by `elastic_fun(x, map_a, map_b)`.
- `bdb_outcome_adj_bin()` and `bdb_outcome_adj_cont()` have been removed.
- **Output restructured.** `mean.theta`/`var.theta` are now
  `post_mean_theta`/`post_var_theta`; `gamma` is replaced by `gamma_base`
  (baseline discounting) and `gamma_base_out` (baseline and outcome
  discounting). New elements are `data_trim`, `N0_trim`, `ps_model`,
  `ps_breaks`, `A_borrow`, `w_s` and `theta_m`. `theta_s` is now an
  `n_strata`-by-3 matrix and `theta_o` a 1-by-3 matrix, both with columns
  `mean`, `lower`, `upper`.
- **Plots** are returned as top-level list elements (`estimate_plot`,
  `ps_group_plot`, `ps_group_stratum_plot`, `density_plot` or
  `density_trace_acf_plot`, `overall_mean_plot`, `stratum_mean_plot`) instead
  of in a `plots` element.

### Methodological changes

- The posterior predictive probability is `P(new <= observed)`, estimated from
  `n_ppp` samples. For binary outcomes the mid-p value (average of the strict
  and non-strict comparison) is used.
- Strata are formed from the quantiles of the current-trial propensity scores
  over their observed range; every stratum must contain at least two subjects
  per arm and data source.
- Stratum-specific posteriors are pooled with weights `w_s = n_s1 / n_1`
  instead of with equal weights.
- `A_eff` is the number of external subjects actually borrowed,
  `sum(n_s0 * gamma_base_out)`.
- Continuous outcomes: the power-prior likelihood includes the normalizing
  constant, and a burn-in (`n_burn`) follows the adaptation phase (`n_adapt`).
- Observed means use t-based confidence intervals; exact (Clopper-Pearson)
  intervals handle zero and full event counts.

### Other changes

- Comprehensive input validation and informative errors and warnings (e.g.
  propensity score model convergence, requested borrowing exceeding the
  available external subjects, large outcome differences under
  `risk_diff = TRUE`).
- Reproducible results through the `seed` argument (R and JAGS). The
  caller's random number generator state is restored on exit.
- The continuous-outcome functions no longer print JAGS progress output.
- `bdb_single_arm_cont()` now also warns when a large outcome difference
  strongly limits borrowing under `risk_diff = TRUE`, and the prior warning of
  the continuous functions takes the scale of the outcome mean into account.
- A warning is issued when the overlap coefficient cannot be computed in a
  stratum (no external data are borrowed there).
- Four simulated data sets are included: `data_bin`, `data_bin_sat`,
  `data_cont` and `data_cont_sat`.

## BayesDDB 0.3.0

This release implements a round of methodological review comments. Two of the
changes alter numerical results; see "Breaking changes" below.

### Breaking changes

- **`gamma_s` now applies the mapped PPP inside `min()`.** The discounting
  factor is now
  `gamma_s = min(1, A * r_s * g(PPP_s) / n_0s)`
  where previously it was `min(1, A * r_s / n_0s) * g(PPP_s)`. Outcome
  discounting is therefore applied to the requested borrowing *before* it is
  capped at the number of external subjects available in the stratum. Results
  are unchanged whenever `A * r_s / n_0s <= 1`; where the requested borrowing
  exceeds the available external subjects, the new form borrows more. The
  `gamma` element of the returned list now carries the outcome discounting,
  so it should no longer be multiplied by `gPPP` downstream.
- **Default `alpha` and `beta` changed from 0.5 to 1 in the binary functions**
  (`bdb_single_arm_bin()`, `bdb_two_arms_bin()`, `bdb_outcome_adj_bin()`).
  `Beta(1, 1)` is the uniform, non-informative prior. The continuous functions
  keep `alpha = beta = 0.01`, since there these are the shape and rate of a
  Gamma prior on the precision, for which `0.01` is the vague choice and
  `Gamma(1, 1)` would be substantially more informative.

### Other changes

- **`a` is now allowed to be any value greater than 0**, rather than being
  restricted to `a >= 1`. Values in `(0, 1)` give a stricter elastic function
  and are mathematically well behaved: `g` remains in `[0, 1]` and still equals
  1 at `PPP = 0.5`.
- **Documentation of `a` corrected.** The previous description had the
  direction backwards. Larger `a` makes the method *more lenient* towards
  outcome differences, not more aggressive. The guidance is now
  `0 < a <= 1` strict, `1 < a <= 5` mild, `a > 10` lenient.
- The vignette gained a numerical table demonstrating the leniency direction,
  and the Step 6 formula was corrected to match the implementation.

## BayesDDB 0.2.0

### New features

- `elastic_function()` gains a power exponent shape parameter `p`, giving the
  two-parameter elastic function
  `g(PPP) = atan(a * sin(pi * PPP))^p / atan(a)^p`. Parameter `a` controls the
  squareness of the curve and `p` controls the sharpness of the peak and the
  steepness of the shoulders. The default `p = 1` reproduces the original
  one-parameter function exactly, so existing results are unchanged.
- `p` is also available on all analysis functions, which now perform the
  outcome-adjustment mapping by calling `elastic_function()` rather than
  repeating the formula in each function.
- `bdb_single_arm_bin()` and `bdb_single_arm_cont()` gain a `risk_diff`
  argument. With `risk_diff = TRUE` the estimand becomes the difference between
  the trial arm and the discounted external data, instead of the trial pooled
  with the discounted external data. The default `risk_diff = FALSE` preserves
  previous behavior.
- New `bdb_outcome_adj_bin()` and `bdb_outcome_adj_cont()` apply the
  outcome-adjustment component alone to a two-arm trial whose control arm is
  augmented with external data. No propensity score model, trimming or
  stratification is performed, so no covariates are required. They return the
  treatment effect in `Estimate` together with `PPP`, `gPPP`, `gamma` and
  `A_eff`.

### Notes

- When trial and external outcomes are strongly incompatible, the elastic
  function drives `gPPP` towards 0 and switches borrowing off. Under
  `risk_diff = TRUE` in `bdb_single_arm_cont()` the external arm is then
  dominated by its diffuse prior and the resulting difference is uninformative.
  Inspect `gPPP` before interpreting the difference.

## BayesDDB 0.1.0

- Initial release.
- Bayesian dynamic borrowing with double adjustment (baseline + outcome).
- Single-arm and two-arm designs; binary and continuous outcomes.
- Propensity score stratification and arctangent elastic function for discounting.
- Functions: `bdb_single_arm_bin`, `bdb_single_arm_cont`, `bdb_two_arms_bin`, `bdb_two_arms_cont`, `elastic_function`.
