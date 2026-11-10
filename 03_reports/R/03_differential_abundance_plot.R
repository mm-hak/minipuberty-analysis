# Plot FDR-significant relative-abundance MaAsLin3 coefficients.
show_significant_relative_da <- function(main_results, interaction_results,
                                         taxonomic_level) {
  results <- dplyr::bind_rows(
    dplyr::mutate(main_results, analysis = "Main-effect model"),
    dplyr::mutate(interaction_results, analysis = "Interaction model")
  ) |>
    dplyr::filter(
      is.finite(coef),
      is.finite(stderr),
      !is.na(qval_joint),
      qval_joint < 0.05,
      !is.na(qval_individual),
      qval_individual < 0.05
    ) |>
    dplyr::distinct(
      feature, source_variable, name, model, coef, stderr,
      qval_individual, qval_joint, N, N_not_zero, analysis,
      .keep_all = TRUE
    ) |>
    dplyr::mutate(
      conf_low = coef - 1.96 * stderr,
      conf_high = coef + 1.96 * stderr,
      result_label = paste0(
        feature, " | ", name, "\n",
        source_variable,
        "; joint q = ", formatC(qval_joint, digits = 2, format = "g"),
        "; component q = ", formatC(qval_individual, digits = 2, format = "g")
      )
    )

  if (nrow(results) == 0) {
    cat(
      "No relative-abundance differential-abundance coefficient met both ",
      "the joint and component FDR thresholds (q < 0.05) at the ",
      taxonomic_level, " level.\n"
    )
    return(invisible(NULL))
  }

  plot <- ggplot2::ggplot(
    results,
    ggplot2::aes(
      x = coef,
      y = stats::reorder(result_label, coef),
      color = analysis,
      shape = model
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = conf_low, xmax = conf_high),
      height = 0.15,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      title = paste0(
        "FDR-significant relative-abundance results: ", taxonomic_level
      ),
      subtitle = paste0(
        "Points are MaAsLin3 coefficients; bars are coefficient ± 1.96 × SE. ",
        "Both joint and component q-values are < 0.05."
      ),
      x = "MaAsLin3 coefficient (95% Wald interval)",
      y = NULL,
      color = "Analysis",
      shape = "Model component"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )

  print(plot)
  invisible(results)
}

