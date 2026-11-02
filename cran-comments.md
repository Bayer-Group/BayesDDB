# CRAN submission comments

## Package

BayesDDB implements Bayesian dynamic borrowing with double adjustment for clinical trials augmented by external data.

## Dependencies

- **Imports:** dplyr, ggplot2, ggpubr, overlapping, stats (required at runtime).
- **Suggests:** rjags — required only for continuous-outcome functions (`bdb_single_arm_cont`, `bdb_two_arms_cont`). The system dependency JAGS (https://mcmc-jags.sourceforge.io/) must be installed separately for those functions. Binary-outcome functions and `elastic_fun` do not require rjags or JAGS. **Check without rjags:** Because rjags/JAGS may be unavailable on some systems (e.g. CRAN build servers), a complete check can be run with `_R_CHECK_FORCE_SUGGESTS_=false` so that missing Suggested packages do not cause an error.

## Examples

- Examples of the four main functions are wrapped in `\donttest{}`. The continuous-outcome examples are additionally guarded by `requireNamespace("rjags", quietly = TRUE)` and take about 15-20 seconds each.

## Testing

- `testthat` tests cover `elastic_fun`, the internal helpers and input checks, the four main `bdb_*` functions, and a regression test against reference results.
- Tests that require rjags/JAGS use `skip_if_not_installed("rjags")` and `skip_on_cran()` so CRAN checks pass without JAGS.

## Test environments

- Ubuntu 22.04.5 LTS, R 4.5.1: `R CMD check --as-cran` with rjags/JAGS 4.3.0 installed, and with rjags unavailable (`_R_CHECK_FORCE_SUGGESTS_=false`) (version 0.4.0).
- Ubuntu 22.04.5 LTS, R 4.4.3: `R CMD check --as-cran` with rjags 4.17 / JAGS 4.3.0 installed (version 0.5.0).

## R CMD check results

0 errors | 0 warnings | 1 note

- NOTE (CRAN incoming feasibility): new submission; the GitHub URLs in `URL` and `BugReports` currently return 404 because the repository is not yet public. They will resolve once the repository is published.

(Apart from messages caused by the local check machine: WARNING that `qpdf` is not installed, NOTE "unable to verify current time" and NOTE that HTML Tidy is not installed.)

## Notes

- No compiler or system libraries required beyond R and the listed R packages.
- Vignette builds with knitr/rmarkdown (Suggests).
