#' Overall observed summaries of external and current data
#'
#' @keywords internal
#' @noRd
plot_group_overall <- function(
  external,
  internal,
  internal_label,
  alpha_level,
  value_label,
  title
) {
  # Overall observed summary for the external and internal data
  #
  # INPUT
  # external: Length-3 vector or 1 by 3 matrix: mean, lower, upper for the external data.
  # internal: Same, for the current trial.
  # internal_label: "SAT" (single arm) or "ICA" (two arms).
  # alpha_level: Significance level, used in the axis label.
  # value_label: "Event rate" (binary) or "Sample Mean" (continuous).
  # title: Plot title.
  # OUTPUT
  # ggplot object.
  ################################################################################
  label <- c("ECA", internal_label)
  mean <- c(external[1], internal[1])
  lower <- c(external[2], internal[2])
  upper <- c(external[3], internal[3])
  group <- c("ECA", internal_label)

  df <- data.frame(label, mean, lower, upper, group)

  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels = rev(df$label))
  df$group <- factor(df$group, levels = rev(df$group))

  fig <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(x = label, y = mean, ymin = lower, ymax = upper, fill = group)
  ) +
    ggplot2::geom_pointrange(linewidth = 1) +
    ggplot2::geom_point(size = 3, shape = 21, stroke = 0.25) +
    ggplot2::coord_flip() + # flip coordinates (puts labels on y axis)
    ggplot2::xlab("") +
    ggplot2::ylab(paste(
      value_label,
      " [",
      100 * (1 - alpha_level),
      "% CI]",
      sep = ""
    )) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() # use a white background
  return(fig)
}
