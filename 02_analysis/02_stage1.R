# Supplementary growth analyses in the blood-sample subset for growth analyses.
# Results from this script are not reported in the manuscript.
# Run via workflow.R or with 02_analysis as the working directory.

source("analysis_helpers.R", local = FALSE)
check_and_load_packages(c("tidyverse", "broom", "rstatix", "patchwork"))

input_path <- file.path("..", "data", "data_reduced.Rds")
output_root <- file.path("output", "stage1")
output_paths <- list(
  tables = file.path(output_root, "tables"),
  figures = file.path(output_root, "figures"),
  models = file.path(output_root, "models"),
  diagnostics = file.path(output_root, "diagnostics"),
  objects = file.path(output_root, "objects")
)
ensure_output_directories(output_paths)
if (!file.exists(input_path)) stop("Processed input not found: ", input_path, call. = FALSE)

data <- readRDS(input_path)

Categoric_variables <- c(
  "Pregnancy.antibiotic", "Labour.antibiotic", "Child.antibiotic.under.7.days",
  "Child.antibiotic.7.days.to.3M", "early_life_ab", "Intervention", "Birth.mode"
)
Continous_pred <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy", "Birth.height", "Birth.weight"
)
outcomes <- c(
  "Growth.velocity.0.3M", "Growth.velocity.3.6M",
  "Weight.gain.0.3M", "Weight.gain.3.6M"
)
outcomes_log <- paste0(outcomes, "_log")
outcomes_list <- list(
  "0.3M" = c("Growth.velocity.0.3M", "Weight.gain.0.3M"),
  "3.6M" = c("Growth.velocity.3.6M", "Weight.gain.3.6M")
)

required <- unique(c(Categoric_variables, Continous_pred, outcomes, outcomes_log,
                     "Growth.velocity.0.6M", "Weight.gain.0.6M"))
missing <- setdiff(required, names(data))
if (length(missing)) stop("Required variable(s) missing: ", paste(missing, collapse = ", "))

# Outcome correlations -------------------------------------------------------
growth_vars <- c("Growth.velocity.0.3M", "Growth.velocity.3.6M", "Growth.velocity.0.6M")
weight_vars <- c("Weight.gain.0.3M", "Weight.gain.3.6M", "Weight.gain.0.6M")
growth_cor <- stats::cor(data[, growth_vars], use = "pairwise.complete.obs")
weight_cor <- stats::cor(data[, weight_vars], use = "pairwise.complete.obs")
save_table(tibble::rownames_to_column(as.data.frame(growth_cor), "variable"),
           file.path(output_paths$tables, "growth_velocity_correlations.csv"))
save_table(tibble::rownames_to_column(as.data.frame(weight_cor), "variable"),
           file.path(output_paths$tables, "weight_gain_correlations.csv"))

# Distribution and assumption plots -----------------------------------------
distribution_plot <- c(outcomes, outcomes_log) |>
  purrr::set_names() |>
  purrr::map(function(variable) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable]])) +
      ggplot2::geom_histogram(fill = "darkblue", bins = 9, color = "white") +
      ggplot2::labs(title = variable, x = NULL, y = "Frequency") +
      ggplot2::theme_minimal()
  }) |>
  patchwork::wrap_plots(ncol = 2)
save_plot(distribution_plot, file.path(output_paths$figures, "outcome_distributions_original_and_log"), 10, 12)

histograms <- list()
scatterplots <- list()
variance_boxplots <- list()
for (predictor in Categoric_variables) {
  histograms[[predictor]] <- plot_histograms_by_group(data, outcomes, predictor)
  variance_boxplots[[predictor]] <- patchwork::wrap_plots(lapply(outcomes, function(outcome) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[predictor]], y = .data[[outcome]])) +
      ggplot2::geom_boxplot(outlier.shape = 21) + ggplot2::theme_minimal() +
      ggplot2::labs(x = predictor, y = outcome)
  }), ncol = 2)
  stem <- safe_file_name(predictor)
  save_plot(histograms[[predictor]], file.path(output_paths$figures, paste0("histograms_by_", stem)), 10, 8)
  save_plot(variance_boxplots[[predictor]], file.path(output_paths$figures, paste0("variance_boxplots_", stem)), 10, 8)
}
for (predictor in Continous_pred) {
  scatterplots[[predictor]] <- patchwork::wrap_plots(
    lapply(outcomes, function(outcome) plot_pair(predictor, outcome, data)), ncol = 2
  )
  save_plot(scatterplots[[predictor]],
            file.path(output_paths$figures, paste0("scatter_", safe_file_name(predictor))), 10, 8)
}

# Mann-Whitney U tests --------------------------------------------------------
wilcoxon <- list()
wilcoxon_plots <- list()
for (period in names(outcomes_list)) {
  wilcoxon[[period]] <- make_wilcox_table(data, outcomes_list[[period]], Categoric_variables)
  save_table(wilcoxon[[period]], file.path(output_paths$tables, paste0("wilcoxon_", period, ".csv")))
  wilcoxon_plots[[period]] <- make_significant_boxplots(wilcoxon[[period]], data)
  save_plot(wilcoxon_plots[[period]], file.path(output_paths$figures, paste0("wilcoxon_selected_", period)), 8, 10)
}

# Continuous-predictor linear models -----------------------------------------
linear_models <- list()
estimate_plots <- list()
for (period in names(outcomes_list)) {
  linear_models[[period]] <- fit_linear_model_grid(data, outcomes_list[[period]], Continous_pred)
  save_table(linear_models[[period]]$results,
             file.path(output_paths$tables, paste0("linear_models_", period, ".csv")))
  saveRDS(linear_models[[period]]$models,
          file.path(output_paths$models, paste0("linear_models_", period, ".rds")))
  save_model_diagnostics(linear_models[[period]]$models,
                         file.path(output_paths$diagnostics, period))
  estimate_plots[[period]] <- plot_model_estimates(linear_models[[period]]$results)
  save_plot(estimate_plots[[period]],
            file.path(output_paths$figures, paste0("linear_model_estimates_", period)), 9, 7)
}

saveRDS(list(correlations = list(growth = growth_cor, weight = weight_cor),
             wilcoxon = wilcoxon, linear_models = linear_models),
        file.path(output_paths$objects, "stage1_results.rds"))
saveRDS(list(distributions = distribution_plot, histograms = histograms,
             scatterplots = scatterplots, variance_boxplots = variance_boxplots,
             wilcoxon = wilcoxon_plots, estimates = estimate_plots),
        file.path(output_paths$objects, "stage1_plots.rds"))

message("Stage 1 analysis complete: ", output_root)
