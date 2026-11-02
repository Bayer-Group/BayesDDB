#' Propensity score density by data source and stratum
#'
#' @keywords internal
#' @noRd
plot_ps_by_stratum <- function(prd, stratum, z) {
  # Propensity score density by data source, faceted by stratum
  #
  # INPUT
  # prd: Propensity scores of the retained subjects.
  # stratum: Stratum index per subject.
  # z: External-data indicator.
  # OUTPUT
  # ggplot object.
  ################################################################################
  df <- data.frame(
    prd = prd,
    stratum = stratum,
    control = ifelse(z == 1, "Current", "External")
  )
  ggplot2::ggplot(df, ggplot2::aes(x = prd, fill = control)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::facet_wrap(~stratum, scales = "free") +
    ggplot2::labs(
      title = "Propensity Score Distribution by Group and Stratum",
      x = "Propensity Score",
      y = "Density",
      fill = "Group"
    ) +
    ggplot2::theme_minimal()
}
