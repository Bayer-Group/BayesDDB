#' Observed summaries of external and current data by stratum
#'
#' @keywords internal
#' @noRd
plot_group_by_stratum <- function(
  external,
  internal,
  internal_label,
  n_strata,
  alpha_level,
  value_label,
  title
) {
  # Observed summary for the external and internal data, by stratum
  #
  # INPUT:
  # external: n_strata-by-3 matrix for the external data.
  # internal: n_strata-by-3 matrix for the current trial.
  # internal_label: "SAT" (single arm) or "ICA" (two arms).
  # n_strata: Number of strata.
  # alpha_level: Significance level, used in the axis label.
  # value_label: "Event rate" (binary) or "Sample Mean" (continuous).
  # title: Plot title.
  # OUTPUT
  # ggplot object.
  ##############################################################################

  label <- rep(stratum_labels(n_strata), 2)
  group <- c(rep("ECA", n_strata), rep(internal_label, n_strata))
  rate <- c(external[, 1], internal[, 1])
  lower <- c(external[, 2], internal[, 2])
  upper <- c(external[, 3], internal[, 3])

  df <- data.frame(group, label, rate, lower, upper)

  # Reverse the factor level ordering for labels after coord_flip()
  df$label <- factor(df$label, levels = rev(df$label)[1:n_strata])
  df$group <- factor(df$group, levels = rev(unique(df$group)))

  axis_text <- paste(
    value_label,
    " [",
    100 * (1 - alpha_level),
    "% CI]",
    sep = ""
  )

  dodge <- ggplot2::position_dodge(width = 0.25)

  fig <- ggplot2::ggplot(
    data = df,
    ggplot2::aes(x = label, y = rate, ymin = lower, ymax = upper, fill = group)
  ) +
    ggplot2::geom_pointrange(linewidth = 1, position = dodge) +
    ggplot2::geom_point(size = 3, shape = 21, stroke = 0.25, position = dodge) +
    ggplot2::scale_x_discrete(name = " ") +
    ggplot2::scale_y_continuous(name = axis_text) +
    ggplot2::coord_flip() +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw()
  return(fig)
}
