#' Stratum axis labels
#'
#' @keywords internal
#' @noRd
stratum_labels <- function(n_strata) {
  # Build the per-stratum axis labels used by the plotting helpers
  # INPUT
  # n_strata: number of strata (positive whole number)
  # OUTPUT
  # character vector of length n_strata: "Stratum 1", "Stratum 2", ...
  paste("Stratum", seq_len(n_strata))
}
