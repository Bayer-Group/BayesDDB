#' Simulated two-arm trial with binary outcome
#'
#' Simulated data set of a randomized trial (2:1, treatment vs control) whose
#' control arm is augmented with external controls, with a binary outcome.
#'
#' @format A data frame with 3300 rows (300 current-trial and 3000 external
#'   subjects) and 13 variables:
#' \describe{
#'   \item{X1, X2, X3, X4}{Binary baseline covariates (0/1).}
#'   \item{X5, X6, X7, X8, X9, X10}{Continuous baseline covariates.}
#'   \item{Z}{Data source: 1 for current trial, 0 for external data.}
#'   \item{TRT}{Treatment: 1 for active treatment, 0 for control.}
#'   \item{Y}{Binary outcome (0/1).}
#' }
#' @source Simulated data.
#' @seealso \code{\link{bdb_two_arms_bin}}
"data_bin"

#' Simulated single-arm trial with binary outcome
#'
#' Simulated data set of a single-arm trial augmented with external data, with a
#' binary outcome.
#'
#' @format A data frame with 2200 rows (200 current-trial and 2000 external
#'   subjects) and 12 variables:
#' \describe{
#'   \item{X1, X2, X3, X4}{Binary baseline covariates (0/1).}
#'   \item{X5, X6, X7, X8, X9, X10}{Continuous baseline covariates.}
#'   \item{Z}{Data source: 1 for current trial, 0 for external data.}
#'   \item{Y}{Binary outcome (0/1).}
#' }
#' @source Simulated data.
#' @seealso \code{\link{bdb_single_arm_bin}}
"data_bin_sat"

#' Simulated two-arm trial with continuous outcome
#'
#' Simulated data set of a randomized trial (2:1, treatment vs control) whose
#' control arm is augmented with external controls, with a continuous outcome.
#'
#' @format A data frame with 3600 rows (600 current-trial and 3000 external
#'   subjects) and 13 variables:
#' \describe{
#'   \item{X1, X2, X3, X4}{Binary baseline covariates (0/1).}
#'   \item{X5, X6, X7, X8, X9, X10}{Continuous baseline covariates.}
#'   \item{Z}{Data source: 1 for current trial, 0 for external data.}
#'   \item{TRT}{Treatment: 1 for active treatment, 0 for control.}
#'   \item{Y}{Continuous outcome.}
#' }
#' @source Simulated data.
#' @seealso \code{\link{bdb_two_arms_cont}}
"data_cont"

#' Simulated single-arm trial with continuous outcome
#'
#' Simulated data set of a single-arm trial augmented with external data, with a
#' continuous outcome.
#'
#' @format A data frame with 3100 rows (100 current-trial and 3000 external
#'   subjects) and 12 variables:
#' \describe{
#'   \item{X1, X2, X3, X4}{Binary baseline covariates (0/1).}
#'   \item{X5, X6, X7, X8, X9, X10}{Continuous baseline covariates.}
#'   \item{Z}{Data source: 1 for current trial, 0 for external data.}
#'   \item{Y}{Continuous outcome.}
#' }
#' @source Simulated data.
#' @seealso \code{\link{bdb_single_arm_cont}}
"data_cont_sat"
