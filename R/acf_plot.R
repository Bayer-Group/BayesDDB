#' Autocorrelation plot of posterior draws
#'
#' @keywords internal
#' @noRd
acf_plot <- function(theta_m, alpha_level) {
  # Autocorrelation function for trace plot
  #
  # INPUT
  # theta_m: vector of posterior draws.
  # alpha_level: significance level
  # OUTPUT
  # ggplot object
  #############################################
  n <- length(theta_m)
  acf_dat <- stats::acf(theta_m, plot = FALSE)
  df <- with(acf_dat, data.frame(lag, acf))
  ciline <- stats::qnorm(1 - alpha_level / 2) / sqrt(n)
  fig <- ggplot2::ggplot(data = df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(x = lag, ymax = acf, ymin = 0),
      width = 0
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = ciline),
      linetype = 2,
      color = "darkblue"
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = -ciline),
      linetype = 2,
      color = "darkblue"
    )
  return(fig)
}
