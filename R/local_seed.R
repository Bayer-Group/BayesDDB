#' Set the RNG seed for the duration of the calling function
#'
#' @keywords internal
#' @noRd
local_seed <- function(seed, envir = parent.frame()) {
  # Seeds the pseudo-random number generator and restores the caller's
  # previous RNG state when the calling function exits.
  # INPUT
  # seed: non-negative whole number, or NULL (no seed is set)
  # envir: frame of the calling function whose exit restores the state
  # OUTPUT
  # invisible NULL
  if (is.null(seed)) {
    return(invisible(NULL))
  }
  old_seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  restore <- function() {
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }
  do.call(on.exit, list(as.call(list(restore)), add = TRUE), envir = envir)
  set.seed(seed)
  invisible(NULL)
}
