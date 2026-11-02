#' Density of the pooled posterior draws
#'
#' @keywords internal
#' @noRd
plot_posterior_density <- function(theta_m) {
  # Density of the pooled posterior draws
  #
  # INPUT
  # theta_m: vector of pooled posterior draws.
  # OUTPUT
  # ggplot object.
  #############################################
  df <- data.frame(x = theta_m)
  ggplot2::ggplot(df, ggplot2::aes(x)) +
    ggplot2::geom_density(alpha = 0.7) +
    ggplot2::labs(
      title = "Density of posterior samples ",
      x = "theta",
      y = "Density"
    )
}
