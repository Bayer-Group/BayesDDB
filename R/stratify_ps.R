#' Stratify subjects by propensity score
#'
#' @keywords internal
#' @noRd
stratify_ps <- function(prd_all, prd_current, n_strata) {
  # INPUT
  # prd_all: Propensity scores of every retained subject.
  # prd_current: Propensity scores of the current-trial subjects, which
  #                define the quantile breaks.
  # n_strata: Number of strata to stratify propensity scores
  # OUTPUT
  # List with stratum (integer vector) and propensity score breaks.
  if (n_strata == 1) {
    # after trimming prd_reference does not span [0,1], only its observed range
    q <- range(prd_current)
    if (q[1] == q[2]) {
      stop(
        "All current-trial propensity scores are identical; stratification is not possible.",
        call. = FALSE
      )
    }
    stratum <- rep(1L, length(prd_all))
  } else {
    q <- stats::quantile(
      prd_current,
      probs = seq(0, 1, length.out = n_strata + 1),
      type = 7,
      names = FALSE
    )
    # Check uniqueness of breaks
    if (anyDuplicated(q)) {
      stop(
        "Propensity-score quantile breaks are not unique; reduce n_strata or revise PS model.",
        call. = FALSE
      )
    }
    stratum <- cut(
      prd_all,
      breaks = q,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    )
  }

  # Check for valid strata
  if (anyNA(stratum)) {
    stop("Failed to assign some observations to strata", call. = FALSE)
  }

  list(stratum = stratum, breaks = q)
}
