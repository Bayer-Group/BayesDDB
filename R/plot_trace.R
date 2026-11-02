#' Trace plot of posterior draws
#'
#' @keywords internal
#' @noRd
plot_trace <- function(theta_m) {
  # Trace plot of posterior draws
  #
  # INPUT
  # theta_m: vector of posterior draws.
  # OUTPUT
  # ggplot object: trace plot
  #############################################
  n_iter <- length(theta_m)
  # Trace plot of mean difference
  df <- data.frame(iteration = 1:n_iter, theta = theta_m)
  fig <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = theta)) +
    ggplot2::geom_line(
      color = "black",
      linewidth = 0.5,
      alpha = 0.9,
      linetype = 1
    ) +
    ggplot2::labs(title = "Trace Plot", x = "Iteration", y = "theta")
  return(fig)
}
