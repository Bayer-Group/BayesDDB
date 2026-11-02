#' Forest plot of the per-stratum and overall estimate
#'
#' @keywords internal
#' @noRd
plot_forest_strata <- function(
  theta_s,
  theta_o,
  n_strata,
  alpha_level,
  title,
  digits,
  hline
) {
  # Forest plot of the per-stratum and overall estimate
  #
  # INPUT
  # theta_s: n_strata-by-3 matrix of mean, lower, upper.
  # theta_o: Overall mean, lower, upper.
  # n_strata: Number of strata.
  # alpha_level: Significance level, used in the axis label.
  # title: Plot title.
  # digits: Rounding for the annotation
  # hline: Draw the dashed zero reference line.
  # OUTPUT
  # ggplot object.
  ################################################################################

  label <- c(stratum_labels(n_strata), "Overall")
  mean <- c(theta_s[, 1], theta_o[1])
  lower <- c(theta_s[, 2], theta_o[2])
  upper <- c(theta_s[, 3], theta_o[3])

  df <- data.frame(label, mean, lower, upper)

  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels = rev(df$label))

  # Dynamic position
  x_pos <- which(levels(df$label) == "Overall") + 0.2

  fig <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(x = label, y = mean, ymin = lower, ymax = upper)
  ) +
    ggplot2::geom_pointrange(linewidth = 1)

  # add a dotted line at x=0 after flip
  if (hline) {
    fig <- fig + ggplot2::geom_hline(yintercept = 0, lty = 2)
  }

  fig +
    ggplot2::coord_flip() + # flip coordinates (puts labels on y axis)
    ggplot2::annotate(
      "text",
      x = x_pos,
      y = theta_o[1],
      label = paste(
        round(theta_o[1], digits),
        " [",
        round(theta_o[2], digits),
        ", ",
        round(theta_o[3], digits),
        "]",
        sep = ""
      )
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab(paste(
      "Posterior Mean [",
      100 * (1 - alpha_level),
      "% CrI]",
      sep = ""
    )) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() # use a white background
}
