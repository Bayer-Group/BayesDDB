#' Propensity score overlapping coefficient by stratum
#'
#' @keywords internal
#' @noRd
overlap_by_stratum <- function(prd, z, stratum, n_strata) {
  # Propensity score overlap coefficient per stratum
  # INPUT
  # prd: Propensity scores of every retained subject.
  # z: External-data indicator (1 current, 0 external).
  # stratum: Integer stratum index per subject.
  # n_strata: Number of strata.
  # OUTPUT
  # Double vector of overlap coefficients, one per stratum.
  v <- rep(0, n_strata)
  failed <- rep(FALSE, n_strata)
  for (i in seq_len(n_strata)) {
    # overlapping probability of stratum i
    temp_lst <- list(
      P0_s = prd[z == 0 & stratum == i],
      P1_s = prd[z == 1 & stratum == i]
    )
    v[i] <- tryCatch(
      as.numeric(overlapping::ovmult(x = temp_lst)),
      error = function(e) NA_real_
    )
    if (!is.finite(v[i])) {
      v[i] <- 0
      failed[i] <- TRUE
    }
  }

  if (!is.finite(sum(v)) || sum(v) <= .Machine$double.eps) {
    stop(
      "All overlap coefficients are zero or non-finite; borrowing weights cannot be normalized",
      call. = FALSE
    )
  }
  if (any(failed)) {
    warning(
      "Overlap coefficient could not be computed in stratum ",
      paste(which(failed), collapse = ", "),
      "; it is set to 0, so no external data are borrowed there",
      call. = FALSE
    )
  }
  v
}
