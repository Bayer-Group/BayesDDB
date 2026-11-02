#' Propensity score density by data source
#'
#' @keywords internal
#' @noRd
plot_ps_by_group <- function(prd, z) {
  # Propensity score density by data source
  #
  # INPUT
  # prd: Propensity scores of the retained subjects.
  # z: External-data indicator.
  # OUTPUT
  # ggplot object.
  ###############################################################################
  df <- data.frame(
    x = prd,
    control = ifelse(z == 1, "Current", "External")
  )
  ggplot2::ggplot(df, ggplot2::aes(x = x, fill = control)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::labs(
      title = "Propensity Score Distribution by Group",
      x = "Propensity Score",
      y = "Density",
      fill = "Group"
    )
}
