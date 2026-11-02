## Not run:
# Simulate example data
library(patchwork)

set.seed(123)
N0 <- 80 # external subjects
N1 <- 40 # current trial (single arm)

# Covariates
X <- rbind(matrix(rnorm(N0 * 2), ncol = 2), matrix(rnorm(N1 * 2), ncol = 2))
colnames(X) <- c("X1", "X2")
Z <- c(rep(0, N0), rep(1, N1))
Y <- c(rnorm(N0, 50, 5), rnorm(N1, 52, 5))

# Run BDB analysis (requires JAGS)
result <- bdb_single_arm_cont(
  Y = Y,
  Z = Z,
  X = X,
  n_target_borrow = 20,
  n_strata = 3,
  map_a = 1,
  map_b = 1,
  n_ppp = 1000,
  alpha_level = 0.05,
  n_iter = 1000,
  n_adapt = 500,
  n_burn = 500,
  gamma_shape = 0.01,
  gamma_rate = 0.01,
  tau_mean = 0.00001,
  doubleadj = TRUE,
  risk_diff = FALSE,
  plots = TRUE,
  seed = 123
)

result$stratum_mean_plot +
  result$overall_mean_plot +
  result$density_trace_acf_plot +
  result$ps_group_stratum_plot +
  result$ps_group_plot +
  result$estimate_plot
