# Stage 1 gender-stratified growth analysis.
# Run via workflow.R or with 02_analysis as the working directory.

packages <- c(
  "tidyverse", "broom", "rstatix", "ARTool", "patchwork"
)
source("analysis_helpers.R", local = FALSE)
check_and_load_packages(packages)

input_path <- file.path("..", "data", "g1_gender.Rds")
output_root <- file.path("output", "stage1_gender")
output_paths <- list(
  tables = file.path(output_root, "tables"),
  figures = file.path(output_root, "figures"),
  models = file.path(output_root, "models"),
  diagnostics = file.path(output_root, "diagnostics"),
  objects = file.path(output_root, "objects")
)
ensure_output_directories(output_paths)

if (!file.exists(input_path)) {
  stop("Processed input not found: ", input_path, call. = FALSE)
}
data <- readRDS(input_path)

Categoric_variables <- c(
  "Pregnancy.antibiotic", "Labour.antibiotic",
  "Child.antibiotic.under.7.days", "Child.antibiotic.7.days.to.3M",
  "early_life_ab", "Intervention", "Birth.mode"
)
Continous_variables <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy",
  "Birth.height", "Birth.weight"
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
Wilcoxon_group_variables <- c("Gender", Categoric_variables)

required_variables <- unique(c(
  Wilcoxon_group_variables, Continous_variables, outcomes, outcomes_log
))
missing_variables <- setdiff(required_variables, names(data))
if (length(missing_variables) > 0) {
  stop("Required variable(s) missing: ", paste(missing_variables, collapse = ", "))
}

# Data distributions ----------------------------------------------------------
distribution_plot <- c(outcomes, outcomes_log) |>
  purrr::set_names() |>
  purrr::map(function(variable) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[variable]])) +
      ggplot2::geom_histogram(fill = "darkblue", bins = 9, color = "white") +
      ggplot2::labs(title = variable, x = NULL, y = "Frequency") +
      ggplot2::theme_minimal()
  }) |>
  patchwork::wrap_plots(ncol = 2)
save_plot(
  distribution_plot,
  file.path(output_paths$figures, "outcome_distributions_original_and_log"),
  width = 10, height = 12
)

class_distribution <- make_class_distribution(data, Categoric_variables)
class_distribution_gender <- make_class_distribution(
  data, Categoric_variables, by_gender = TRUE
)
save_table(class_distribution, file.path(output_paths$tables, "class_distribution.csv"))
save_table(
  class_distribution_gender,
  file.path(output_paths$tables, "class_distribution_by_gender.csv")
)

hist_per_group <- list()
hist_per_group_gender <- list()
for (variable in Categoric_variables) {
  hist_per_group[[variable]] <- plot_histograms_by_group(data, outcomes, variable)
  hist_per_group_gender[[variable]] <- plot_histograms_by_group_with_gender(
    data, outcomes, variable
  )
  stem <- safe_file_name(variable)
  save_plot(
    hist_per_group[[variable]],
    file.path(output_paths$figures, paste0("histograms_by_", stem)),
    width = 10, height = 8
  )
  save_plot(
    hist_per_group_gender[[variable]],
    file.path(output_paths$figures, paste0("histograms_by_", stem, "_and_gender")),
    width = 10, height = 8
  )
}

scatter_original <- list()
scatter_original_gender <- list()
scatter_log <- list()
scatter_log_gender <- list()
for (predictor in Continous_variables) {
  scatter_original[[predictor]] <- patchwork::wrap_plots(
    lapply(outcomes, function(outcome) plot_pair(predictor, outcome, data)), ncol = 2
  )
  scatter_original_gender[[predictor]] <- patchwork::wrap_plots(
    lapply(outcomes, function(outcome) plot_pair(predictor, outcome, data, TRUE)),
    ncol = 2, guides = "collect"
  )
  scatter_log[[predictor]] <- patchwork::wrap_plots(
    lapply(outcomes_log, function(outcome) plot_pair(predictor, outcome, data)), ncol = 2
  )
  scatter_log_gender[[predictor]] <- patchwork::wrap_plots(
    lapply(outcomes_log, function(outcome) plot_pair(predictor, outcome, data, TRUE)),
    ncol = 2, guides = "collect"
  )
  stem <- safe_file_name(predictor)
  save_plot(scatter_original[[predictor]], file.path(output_paths$figures, paste0("scatter_", stem)), 10, 8)
  save_plot(scatter_original_gender[[predictor]], file.path(output_paths$figures, paste0("scatter_", stem, "_gender")), 10, 8)
  save_plot(scatter_log[[predictor]], file.path(output_paths$figures, paste0("scatter_", stem, "_log_outcomes")), 10, 8)
  save_plot(scatter_log_gender[[predictor]], file.path(output_paths$figures, paste0("scatter_", stem, "_log_outcomes_gender")), 10, 8)
}

boxplots <- list()
boxplots_grouped <- list()
for (predictor in Categoric_variables) {
  boxplots[[predictor]] <- patchwork::wrap_plots(lapply(outcomes, function(outcome) {
    ggplot2::ggplot(data, ggplot2::aes(x = .data[[predictor]], y = .data[[outcome]])) +
      ggplot2::geom_boxplot(outlier.shape = 21) + ggplot2::theme_minimal() +
      ggplot2::labs(x = predictor, y = outcome)
  }), ncol = 2)
  boxplots_grouped[[predictor]] <- patchwork::wrap_plots(lapply(outcomes, function(outcome) {
    ggplot2::ggplot(
      data, ggplot2::aes(x = .data[[predictor]], y = .data[[outcome]], fill = Gender)
    ) + ggplot2::geom_boxplot(position = ggplot2::position_dodge2(width = 0.8)) +
      ggplot2::theme_minimal() + ggplot2::labs(x = predictor, y = outcome)
  }), ncol = 2)
  stem <- safe_file_name(predictor)
  save_plot(boxplots[[predictor]], file.path(output_paths$figures, paste0("variance_boxplots_", stem)), 10, 8)
  save_plot(boxplots_grouped[[predictor]], file.path(output_paths$figures, paste0("variance_boxplots_", stem, "_gender")), 10, 8)
}

# Mann-Whitney U tests ---------------------------------------------------------
res_tab_Wilcox <- list()
wilcox_plots <- list()
for (time_period in names(outcomes_list)) {
  result <- make_wilcox_table(
    data, outcomes_list[[time_period]], Wilcoxon_group_variables
  )
  res_tab_Wilcox[[time_period]] <- result
  save_table(
    result,
    file.path(output_paths$tables, paste0("wilcoxon_", time_period, ".csv"))
  )
  wilcox_plots[[time_period]] <- make_significant_boxplots(result, data)
  save_plot(
    wilcox_plots[[time_period]],
    file.path(output_paths$figures, paste0("wilcoxon_selected_", time_period)),
    width = 8, height = 10
  )
}

# ART-ANOVA --------------------------------------------------------------------
predictors2 <- c(
  "Intervention", "Birth.mode", "Pregnancy.antibiotic",
  "Labour.antibiotic", "Child.antibiotic.under.7.days",
  "Child.antibiotic.7.days.to.3M", "early_life_ab"
)
art_results_all <- list()
art_failures <- list()
art_plots <- list()
for (time_period in names(outcomes_list)) {
  art <- run_art_models(data, outcomes_list[[time_period]], predictors2, time_period)
  art_results_all[[time_period]] <- art$results
  art_failures[[time_period]] <- art$failures
  save_table(art$results, file.path(output_paths$tables, paste0("art_anova_", time_period, ".csv")))
  if (nrow(art$failures) > 0) {
    save_table(art$failures, file.path(output_paths$tables, paste0("art_anova_failures_", time_period, ".csv")))
  }

  selected <- art$results[art$results$raw_p < 0.05, , drop = FALSE]
  if (nrow(selected) > 0) {
    selected$label <- paste(selected$Outcome, selected$Effect, sep = " | ")
    art_plots[[time_period]] <- ggplot2::ggplot(
      selected, ggplot2::aes(x = effect_size, y = reorder(label, effect_size))
    ) + ggplot2::geom_point(color = "darkblue") + ggplot2::theme_minimal() +
      ggplot2::labs(x = "Partial effect size", y = NULL, title = paste("ART-ANOVA", time_period))
    save_plot(art_plots[[time_period]], file.path(output_paths$figures, paste0("art_anova_selected_", time_period)), 9, 7)
  }
}

# Linear regression ------------------------------------------------------------
results_lm <- list()
results_gender_lm <- list()
lm_estimate_plots <- list()
gender_estimate_plots <- list()
for (time_period in names(outcomes_list)) {
  simple <- fit_linear_model_grid(
    data, outcomes_list[[time_period]], Continous_variables
  )
  gender <- fit_linear_model_grid(
    data, outcomes_list[[time_period]], Continous_variables,
    gender_interaction = TRUE
  )
  results_lm[[time_period]] <- simple
  results_gender_lm[[time_period]] <- gender

  save_table(simple$results, file.path(output_paths$tables, paste0("linear_models_", time_period, ".csv")))
  save_table(gender$results, file.path(output_paths$tables, paste0("linear_models_gender_", time_period, ".csv")))
  saveRDS(simple$models, file.path(output_paths$models, paste0("linear_models_", time_period, ".rds")))
  saveRDS(gender$models, file.path(output_paths$models, paste0("linear_models_gender_", time_period, ".rds")))

  save_model_diagnostics(simple$models, file.path(output_paths$diagnostics, "linear", time_period))
  save_model_diagnostics(gender$models, file.path(output_paths$diagnostics, "gender_interaction", time_period))

  lm_estimate_plots[[time_period]] <- plot_model_estimates(simple$results)
  gender_estimate_plots[[time_period]] <- plot_model_estimates(
    gender$results, interaction_only = TRUE
  )
  save_plot(lm_estimate_plots[[time_period]], file.path(output_paths$figures, paste0("linear_model_estimates_", time_period)), 9, 7)
  save_plot(gender_estimate_plots[[time_period]], file.path(output_paths$figures, paste0("gender_interaction_estimates_", time_period)), 9, 7)
}

# Save complete R objects for validation and report reuse ----------------------
analysis_results <- list(
  variables = list(
    Categoric_variables = Categoric_variables,
    Wilcoxon_group_variables = Wilcoxon_group_variables,
    Continous_variables = Continous_variables,
    outcomes = outcomes,
    outcomes_log = outcomes_log,
    outcomes_list = outcomes_list,
    predictors2 = predictors2
  ),
  class_distribution = class_distribution,
  class_distribution_gender = class_distribution_gender,
  wilcoxon = res_tab_Wilcox,
  art_anova = art_results_all,
  art_failures = art_failures,
  linear_models = results_lm,
  gender_linear_models = results_gender_lm
)
saveRDS(analysis_results, file.path(output_paths$objects, "stage1_gender_results.rds"))
saveRDS(
  list(
    distributions = distribution_plot,
    histograms = hist_per_group,
    histograms_gender = hist_per_group_gender,
    scatter = scatter_original,
    scatter_gender = scatter_original_gender,
    scatter_log = scatter_log,
    scatter_log_gender = scatter_log_gender,
    variance_boxplots = boxplots,
    variance_boxplots_gender = boxplots_grouped,
    wilcoxon = wilcox_plots,
    art = art_plots,
    linear_estimates = lm_estimate_plots,
    gender_estimates = gender_estimate_plots
  ),
  file.path(output_paths$objects, "stage1_gender_plots.rds")
)

message("Stage 1 gender analysis complete: ", output_root)
