#' Trim subjects to the current-trial propensity score range
#'
#' @keywords internal
#' @noRd
trimming <- function(data, prd, prd_current, sort_by_z = FALSE) {
  # INPUT
  # data: data frame with PS scores and Z indicator (1: current, 0: external)
  # prd: propensity scores
  # prd_current: Propensity scores of the current-trial subjects, whose range
  #              defines the retained interval
  # sort_by_z: sort results by Z descending
  # OUTPUT
  # trimmed data frame
  data_trim <- data[dplyr::between(prd, min(prd_current), max(prd_current)), ]
  if (sort_by_z) {
    data_trim <- data_trim[order(data_trim$Z, decreasing = TRUE), ]
  }
  data_trim
}
