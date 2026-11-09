# Descriptive summaries of the processed datasets.
# Run via workflow.R or source this script with 02_analysis as the working directory.

# Load helper functions and packages.
source("analysis_helpers.R", local = FALSE)
check_and_load_packages(c(
  "tidyverse", "SummarizedExperiment", "SingleCellExperiment"
))

# Load processed datasets.
gender_data_path <- file.path(
  "..", "data", "g1_gender.Rds"
)
tse_path <- file.path(
  "..", "data", "tse.Rds"
)
if (!file.exists(gender_data_path) || !file.exists(tse_path)) {
  stop("Processed datasets are missing. Run 01_data_processing/Data.R first.")
}
gender_data <- readRDS(gender_data_path)
tse <- readRDS(tse_path)
tse_metadata <- as.data.frame(SummarizedExperiment::colData(tse))

# Match gender to microbiome samples without changing the TSE.
sample_match <- match(rownames(tse_metadata), rownames(gender_data))
if (anyNA(sample_match)) {
  stop("Some microbiome sample IDs are missing from g1_gender.Rds.")
}
tse_metadata$Gender <- gender_data$Gender[sample_match]

# Define output folders.
output_root <- file.path("output", "data_description")
output_paths <- list(
  tables = file.path(output_root, "tables"),
  figures = file.path(output_root, "figures"),
  objects = file.path(output_root, "objects")
)
ensure_output_directories(output_paths)

# Define descriptive variables.
categorical_variables <- c(
  "Gender", "Pregnancy.antibiotic", "Labour.antibiotic",
  "Child.antibiotic.under.7.days", "Child.antibiotic.7.days.to.3M",
  "early_life_ab", "Intervention", "Birth.mode"
)
growth_variables <- c(
  "Growth.velocity.0.3M", "Growth.velocity.3.6M",
  "Growth.velocity.0.6M", "Weight.gain.0.3M",
  "Weight.gain.3.6M", "Weight.gain.0.6M"
)
background_variables <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy",
  "Birth.height", "Birth.weight"
)
hormone_variables <- c(
  "S.TestoMS", "P.LH", "P.FSH", "P.SHBG",
  "S.TestoMS_log", "P.LH_log", "P.FSH_log", "P.SHBG_log"
)
diversity_variables <- c(
  "dbp_dominance", "faith_diversity", "observed_richness",
  "shannon_diversity"
)

# Summarise categorical variables.
categorical_summary <- function(data, dataset, variables) {
  variables <- intersect(variables, names(data))
  data |>
    dplyr::select(dplyr::all_of(variables)) |>
    tidyr::pivot_longer(
      dplyr::everything(), names_to = "variable", values_to = "level"
    ) |>
    dplyr::mutate(
      level = dplyr::if_else(is.na(level), "Missing", as.character(level))
    ) |>
    dplyr::count(variable, level, name = "n") |>
    dplyr::group_by(variable) |>
    dplyr::mutate(percent = 100 * n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::mutate(dataset = dataset, .before = 1)
}

# Summarise numeric variables.
numeric_summary <- function(data, dataset, variables) {
  variables <- intersect(variables, names(data))
  purrr::map_dfr(variables, function(variable) {
    values <- data[[variable]]
    observed <- values[is.finite(values)]
    data.frame(
      dataset = dataset,
      variable = variable,
      n = length(observed),
      missing = sum(!is.finite(values)),
      mean = mean(observed),
      sd = stats::sd(observed),
      median = stats::median(observed),
      q1 = unname(stats::quantile(observed, 0.25)),
      q3 = unname(stats::quantile(observed, 0.75)),
      minimum = min(observed),
      maximum = max(observed)
    )
  })
}

# Summarise missing values.
missing_summary <- function(data, dataset) {
  tibble::tibble(
    dataset = dataset,
    variable = names(data),
    n = nrow(data),
    missing = vapply(data, function(x) sum(is.na(x)), integer(1)),
    percent_missing = 100 * missing / n
  )
}

# Calculate descriptive correlations.
correlation_summary <- function(data, dataset, variables) {
  variables <- intersect(variables, names(data))
  pairs <- utils::combn(variables, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(pair) {
    complete <- stats::complete.cases(data[, pair])
    data.frame(
      dataset = dataset,
      variable1 = pair[[1]],
      variable2 = pair[[2]],
      n = sum(complete),
      pearson_correlation = stats::cor(
        data[[pair[[1]]]][complete], data[[pair[[2]]]][complete],
        method = "pearson"
      ),
      spearman_correlation = stats::cor(
        data[[pair[[1]]]][complete], data[[pair[[2]]]][complete],
        method = "spearman"
      )
    )
  })
}

# Dataset overview.
dataset_overview <- tibble::tibble(
  dataset = c("Gender and growth dataset", "Microbiome TSE dataset"),
  samples = c(nrow(gender_data), ncol(tse)),
  variables = c(
    ncol(gender_data), ncol(SummarizedExperiment::colData(tse))
  ),
  matched_to_gender_dataset = c(nrow(gender_data), sum(!is.na(sample_match)))
)

taxonomic_levels <- SingleCellExperiment::altExpNames(tse)
taxa_summary <- purrr::map_dfr(taxonomic_levels, function(level) {
  experiment <- SingleCellExperiment::altExp(tse, level)
  tibble::tibble(
    taxonomic_level = level,
    taxa = nrow(experiment),
    samples = ncol(experiment)
  )
})

# Group distributions.
group_distributions <- dplyr::bind_rows(
  categorical_summary(
    gender_data, "Gender and growth dataset", categorical_variables
  ),
  categorical_summary(
    tse_metadata, "Microbiome TSE dataset", categorical_variables
  )
)

# Numeric summaries.
numeric_summaries <- dplyr::bind_rows(
  numeric_summary(
    gender_data, "Gender and growth dataset",
    c(background_variables, growth_variables)
  ),
  numeric_summary(
    tse_metadata, "Microbiome TSE dataset",
    c(background_variables, growth_variables, hormone_variables,
      diversity_variables)
  )
)

# Missing-value summaries.
missing_values <- dplyr::bind_rows(
  missing_summary(gender_data, "Gender and growth dataset"),
  missing_summary(tse_metadata, "Microbiome TSE dataset")
)

# Correlation summaries.
growth_correlations <- dplyr::bind_rows(
  correlation_summary(
    gender_data, "Gender and growth dataset", growth_variables
  ),
  correlation_summary(
    tse_metadata, "Microbiome TSE dataset", growth_variables
  )
)
hormone_correlations <- correlation_summary(
  tse_metadata, "Microbiome TSE dataset", hormone_variables
)

# Gender distribution figure.
gender_plot_data <- group_distributions |>
  dplyr::filter(variable == "Gender", level != "Missing")
gender_distribution_plot <- ggplot2::ggplot(
  gender_plot_data,
  ggplot2::aes(x = level, y = n, fill = dataset)
) +
  ggplot2::geom_col(position = "dodge") +
  ggplot2::geom_text(
    ggplot2::aes(label = paste0(n, " (", round(percent, 1), "%)")),
    position = ggplot2::position_dodge(width = 0.9), vjust = -0.3
  ) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Gender", y = "Number of participants", fill = NULL)
save_plot(
  gender_distribution_plot,
  file.path(output_paths$figures, "gender_distribution"),
  8, 6
)

# Distribution figures.
plot_distributions <- function(data, variables, title) {
  variables <- intersect(variables, names(data))
  data |>
    dplyr::select(dplyr::all_of(variables)) |>
    tidyr::pivot_longer(
      dplyr::everything(), names_to = "variable", values_to = "value"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(bins = 15, fill = "darkblue", color = "white") +
    ggplot2::facet_wrap(~variable, scales = "free", ncol = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = NULL, y = "Count")
}

save_plot(
  plot_distributions(
    gender_data, growth_variables,
    "Growth-variable distributions: gender and growth dataset"
  ),
  file.path(output_paths$figures, "growth_distributions_gender_dataset"),
  10, 10
)
save_plot(
  plot_distributions(
    tse_metadata, growth_variables,
    "Growth-variable distributions: microbiome dataset"
  ),
  file.path(output_paths$figures, "growth_distributions_tse_dataset"),
  10, 10
)
save_plot(
  plot_distributions(
    tse_metadata, hormone_variables[1:4],
    "Hormone distributions: microbiome dataset"
  ),
  file.path(output_paths$figures, "hormone_distributions_tse_dataset"),
  10, 8
)
save_plot(
  plot_distributions(
    tse_metadata, diversity_variables,
    "Diversity distributions: microbiome dataset"
  ),
  file.path(output_paths$figures, "diversity_distributions_tse_dataset"),
  10, 8
)

# Save descriptive tables.
save_table(dataset_overview, file.path(output_paths$tables, "dataset_overview.csv"))
save_table(taxa_summary, file.path(output_paths$tables, "taxa_summary.csv"))
save_table(group_distributions, file.path(output_paths$tables, "group_distributions.csv"))
save_table(numeric_summaries, file.path(output_paths$tables, "numeric_summaries.csv"))
save_table(missing_values, file.path(output_paths$tables, "missing_values.csv"))
save_table(growth_correlations, file.path(output_paths$tables, "growth_correlations.csv"))
save_table(hormone_correlations, file.path(output_paths$tables, "hormone_correlations.csv"))

# Save complete descriptive results.
data_description_results <- list(
  dataset_overview = dataset_overview,
  taxa_summary = taxa_summary,
  group_distributions = group_distributions,
  numeric_summaries = numeric_summaries,
  missing_values = missing_values,
  growth_correlations = growth_correlations,
  hormone_correlations = hormone_correlations
)
saveRDS(
  data_description_results,
  file.path(output_paths$objects, "data_description_results.rds")
)

message("Data description complete: ", output_root)
