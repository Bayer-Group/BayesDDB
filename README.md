# BayesDDB

Bayesian Dynamic Borrowing with Double Adjustment for augmenting clinical trial data with external data.

## Installation

```r
# From CRAN (when available)
install.packages("BayesDDB")

# Development version from GitHub
pak::pkg_install("bayer-group/BayesDDB")
```

## Usage

The package provides four analysis functions and a utility.

Double adjustment (propensity score baseline adjustment + outcome adjustment):

- `bdb_single_arm_bin()` — single-arm trial, binary outcome
- `bdb_single_arm_cont()` — single-arm trial, continuous outcome (requires rjags and JAGS)
- `bdb_two_arms_bin()` — two-arm trial, binary outcome
- `bdb_two_arms_cont()` — two-arm trial, continuous outcome (requires rjags and JAGS)

In the single-arm functions, `risk_diff = TRUE` estimates the difference between
the trial arm and the external data as comparator. Outcome adjustment is not
valid for this estimand, so only baseline adjustment is applied.

Utility:

- `elastic_fun()` — arctangent elastic function for outcome discounting

Minimal example with the elastic function:

```r
library(BayesDDB)
elastic_fun(c(0.1, 0.5, 0.9), map_a = 1, map_b = 1)

# The power parameter map_b sharpens the peak (map_b = 1 is the original form)
elastic_fun(c(0.1, 0.5, 0.9), map_a = 1, map_b = 3)
```

The elastic function has two shape parameters: `map_a` controls the squareness of the
curve and `map_b` controls the sharpness of the peak. Both are available as arguments
on every analysis function.

For full methodology and runnable examples, see the vignette:

```r
vignette("introduction", package = "BayesDDB")
```

To install from source with the vignette built: 
`remotes::install_github("bayer-group/BayesDDB", build_vignettes = TRUE)`.

## License

GPL (>= 3).
