#' Count subjects per stratum
#'
#' @keywords internal
#' @noRd
stratum_counts <- function(stratum, keep, n_strata) {
  # stratum_counts() counts subjects per stratum within a subgroup
  # INPUT
  # stratum: Integer stratum index per subject.
  # keep: Logical vector selecting the subgroup.
  # n_strata: Number of strata.
  # OUTPUT
  # Vector of length n_strata.
  return(tabulate(stratum[keep], nbins = n_strata))
}
