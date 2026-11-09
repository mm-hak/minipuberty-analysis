# Stage 2 analysis of hormone levels and their associations with growth.
# Run via workflow.R or with 02_analysis as the working directory.

source("analysis_helpers.R", local = FALSE)
check_and_load_packages(c("tidyverse", "broom", "rstatix", "patchwork"))

# Load the processed data -----------------------------------------------------
input_path <- file.path("..", "data", "data_reduced.Rds")
output_root <- file.path("output", "stage2")
output_paths <- list(tables = file.path(output_root, "tables"),
                     figures = file.path(output_root, "figures"),
                     models = file.path(output_root, "models"),
                     diagnostics = file.path(output_root, "diagnostics"),
                     objects = file.path(output_root, "objects"))
ensure_output_directories(output_paths)
if (!file.exists(input_path)) stop("Processed input not found: ", input_path, call. = FALSE)
data <- readRDS(input_path)

# Define the analysis variables ----------------------------------------------
Categoric_variables <- c("Pregnancy.antibiotic", "Labour.antibiotic",
  "Child.antibiotic.under.7.days", "Child.antibiotic.7.days.to.3M",
  "early_life_ab", "Intervention", "Birth.mode")
Continous_pred <- c("Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy", "Birth.height", "Birth.weight",
  "Growth.velocity.0.3M", "Growth.velocity.3.6M", "Weight.gain.0.3M", "Weight.gain.3.6M")
outcomes <- c("Growth.velocity.0.3M", "Growth.velocity.3.6M", "Weight.gain.0.3M", "Weight.gain.3.6M")
hormones <- c("S.TestoMS", "P.LH", "P.FSH", "P.SHBG")
hormones_log <- paste0(hormones, "_log")
outcomes_list <- list("0.3M" = c("Growth.velocity.0.3M", "Weight.gain.0.3M"),
                      "3.6M" = c("Growth.velocity.3.6M", "Weight.gain.3.6M"))
missing <- setdiff(unique(c(Categoric_variables, Continous_pred, outcomes, hormones, hormones_log)), names(data))
if (length(missing)) stop("Required variable(s) missing: ", paste(missing, collapse = ", "))

# Hormone distributions -------------------------------------------------------
make_distribution <- function(vars) {
  data |> dplyr::select(dplyr::all_of(vars)) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "variable", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(fill = "darkblue", bins = 9, color = "white") +
    ggplot2::facet_wrap(~variable, scales = "free") + ggplot2::theme_minimal()
}
distribution_original <- make_distribution(hormones)
distribution_log <- make_distribution(hormones_log)
save_plot(distribution_original, file.path(output_paths$figures, "hormone_distributions_original"), 10, 7)
save_plot(distribution_log, file.path(output_paths$figures, "hormone_distributions_log"), 10, 7)

# Exploratory plots -----------------------------------------------------------
histograms <- list(); scatter_predictors <- list(); scatter_growth <- list(); variance_boxplots <- list()
for (predictor in Categoric_variables) {
  histograms[[predictor]] <- list(original = plot_histograms_by_group(data, hormones, predictor),
                                 log = plot_histograms_by_group(data, hormones_log, predictor))
  variance_boxplots[[predictor]] <- patchwork::wrap_plots(lapply(hormones_log, function(hormone) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[predictor]], y = .data[[hormone]])) +
      ggplot2::geom_boxplot(outlier.shape = 21) + ggplot2::theme_minimal()
  }), ncol = 2)
  stem <- safe_file_name(predictor)
  save_plot(histograms[[predictor]]$original, file.path(output_paths$figures, paste0("histograms_original_by_", stem)), 10, 8)
  save_plot(histograms[[predictor]]$log, file.path(output_paths$figures, paste0("histograms_log_by_", stem)), 10, 8)
  save_plot(variance_boxplots[[predictor]], file.path(output_paths$figures, paste0("variance_boxplots_", stem)), 10, 8)
}
for (predictor in Continous_pred) {
  scatter_predictors[[predictor]] <- list(
    original = patchwork::wrap_plots(lapply(hormones, function(h) plot_pair(predictor, h, data)), ncol = 2),
    log = patchwork::wrap_plots(lapply(hormones_log, function(h) plot_pair(predictor, h, data)), ncol = 2))
  stem <- safe_file_name(predictor)
  save_plot(scatter_predictors[[predictor]]$original, file.path(output_paths$figures, paste0("predictor_hormone_", stem)), 10, 8)
  save_plot(scatter_predictors[[predictor]]$log, file.path(output_paths$figures, paste0("predictor_log_hormone_", stem)), 10, 8)
}
for (hormone in hormones_log) {
  scatter_growth[[hormone]] <- patchwork::wrap_plots(lapply(outcomes, function(y) plot_pair(hormone, y, data)), ncol = 2)
  save_plot(scatter_growth[[hormone]], file.path(output_paths$figures, paste0("growth_by_", safe_file_name(hormone))), 10, 8)
}

# Mann-Whitney U tests --------------------------------------------------------
wilcoxon <- make_wilcox_table(data, hormones_log, Categoric_variables)
save_table(wilcoxon, file.path(output_paths$tables, "wilcoxon_hormones.csv"))
wilcoxon_plot <- make_significant_boxplots(wilcoxon, data)
save_plot(wilcoxon_plot, file.path(output_paths$figures, "wilcoxon_hormones_selected"), 8, 10)

# Hormones versus continuous predictors --------------------------------------
hormones_vs_predictors <- fit_linear_model_grid(data, hormones_log, Continous_pred)
save_table(hormones_vs_predictors$results, file.path(output_paths$tables, "linear_hormones_vs_predictors.csv"))
saveRDS(hormones_vs_predictors$models, file.path(output_paths$models, "linear_hormones_vs_predictors.rds"))
save_model_diagnostics(hormones_vs_predictors$models, file.path(output_paths$diagnostics, "hormones_vs_predictors"))
hp_plot <- plot_model_estimates(hormones_vs_predictors$results)
save_plot(hp_plot, file.path(output_paths$figures, "linear_hormones_vs_predictors_estimates"), 9, 7)

# Growth outcomes versus hormones --------------------------------------------
growth_vs_hormones <- list(); gh_plots <- list()
for (period in names(outcomes_list)) {
  growth_vs_hormones[[period]] <- fit_linear_model_grid(data, outcomes_list[[period]], hormones_log)
  save_table(growth_vs_hormones[[period]]$results,
             file.path(output_paths$tables, paste0("linear_growth_vs_hormones_", period, ".csv")))
  saveRDS(growth_vs_hormones[[period]]$models,
          file.path(output_paths$models, paste0("linear_growth_vs_hormones_", period, ".rds")))
  save_model_diagnostics(growth_vs_hormones[[period]]$models,
                         file.path(output_paths$diagnostics, "growth_vs_hormones", period))
  gh_plots[[period]] <- plot_model_estimates(growth_vs_hormones[[period]]$results)
  save_plot(gh_plots[[period]], file.path(output_paths$figures, paste0("linear_growth_vs_hormones_", period)), 9, 7)
}

# Testosterone interaction models -------------------------------------------
interaction_models <- list(
  bmi = stats::lm(S.TestoMS_log ~ Pregnancy.antibiotic * Mother.BMI.before.pregnancy, data = data),
  maternal_age = stats::lm(S.TestoMS_log ~ Pregnancy.antibiotic * Mothers.age, data = data)
)
interaction_table <- purrr::imap_dfr(interaction_models, function(model, name) {
  broom::tidy(model) |> dplyr::mutate(model = name)
}) |> dplyr::filter(term != "(Intercept)") |>
  dplyr::mutate(p_fdr_all = stats::p.adjust(p.value, method = "fdr"))
save_table(interaction_table, file.path(output_paths$tables, "testosterone_interaction_models.csv"))
saveRDS(interaction_models, file.path(output_paths$models, "testosterone_interaction_models.rds"))

interaction_bmi_plot <- ggplot2::ggplot(data, ggplot2::aes(
  x = Mother.BMI.before.pregnancy, y = S.TestoMS_log, color = factor(Pregnancy.antibiotic))) +
  ggplot2::geom_point() + ggplot2::geom_smooth(method = "lm", se = TRUE) + ggplot2::theme_minimal()
interaction_age_plot <- ggplot2::ggplot(data, ggplot2::aes(
  x = Mothers.age, y = S.TestoMS_log, color = factor(Pregnancy.antibiotic))) +
  ggplot2::geom_point() + ggplot2::geom_smooth(method = "lm", se = TRUE) + ggplot2::theme_minimal()
save_plot(interaction_bmi_plot, file.path(output_paths$figures, "testosterone_pregnancy_antibiotic_bmi"), 8, 6)
save_plot(interaction_age_plot, file.path(output_paths$figures, "testosterone_pregnancy_antibiotic_maternal_age"), 8, 6)

# Save the complete Stage 2 results ------------------------------------------
saveRDS(list(wilcoxon = wilcoxon, hormones_vs_predictors = hormones_vs_predictors,
             growth_vs_hormones = growth_vs_hormones, interactions = interaction_table),
        file.path(output_paths$objects, "stage2_results.rds"))
message("Stage 2 analysis complete: ", output_root)
