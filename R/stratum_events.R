#' Count events per stratum
#'
#' @keywords internal
#' @noRd
stratum_events <- function(values, stratum, keep, n_strata) {
  # Sum binary events (0/1) per stratum within a subgroup
  #
  # INPUT
  # values: Numeric/binary vector containing 0s and 1s.
  # stratum: Integer stratum index per subject.
  # keep: Logical vector selecting the subgroup.
  # n_strata: Number of strata.
  #
  # OUTPUT
  # Returns: Numeric vector of length n_strata containing event counts.
  stopifnot(all(values %in% c(0, 1))) # input check
  tabulate(
    bin = stratum[keep & values == 1],
    nbins = n_strata
  )
}
