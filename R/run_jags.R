#' Adapt, burn in and sample from a JAGS model
#'
#' @keywords internal
#' @noRd
run_jags <- function(
  model_string,
  data,
  inits,
  params,
  n_adapt,
  n_burn,
  n_iter,
  seed = NULL
) {
  # Adapt, burn in and sample from a JAGS model
  #
  # INPUT
  # model_string: Model definition as a character scalar.
  # data: Named list handed to JAGS as data.
  # inits: Initial values: a named list or a function returning one.
  # params: Character vector of parameters to monitor.
  # n_adapt: Adaptation iterations.
  # n_burn: Burn-in iterations discarded via update().
  # n_iter: Monitored iterations.
  # seed: seed for pseudo-random generator
  # OUTPUT
  # Data frame of posterior draws, one column per monitored parameter.
  ###############################################################################

  if (!is.null(seed)) {
    if (is.null(inits)) {
      inits <- list()
    }
    inits$.RNG.name <- "base::Mersenne-Twister"
    # Callers add offsets to the user seed, which can exceed the integer range
    inits$.RNG.seed <- as.integer(seed %% .Machine$integer.max)
  }

  con <- textConnection(model_string)
  on.exit(close(con), add = TRUE)

  fit <- rjags::jags.model(
    con,
    data = data,
    inits = inits,
    n.chains = 1,
    n.adapt = n_adapt,
    quiet = TRUE
  )
  stats::update(fit, n_burn, progress.bar = "none")
  draws <- rjags::coda.samples(
    fit,
    variable.names = params,
    n.iter = n_iter,
    thin = 1,
    progress.bar = "none"
  )
  do.call(rbind.data.frame, draws)
}
