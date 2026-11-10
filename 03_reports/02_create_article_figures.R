# Create the final multi-panel figures for the article.
# Run via workflow.R or with 03_reports as the working directory.

# Check required packages.
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "patchwork",
  "SummarizedExperiment", "SingleCellExperiment"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Required package(s) are not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# Define paths.
data_directory <- file.path("..", "data")
analysis_directory <- file.path("..", "02_analysis", "output")
output_directory <- file.path("output", "article_figures")

gender_data_path <- file.path(data_directory, "g1_gender.Rds")
reduced_data_path <- file.path(data_directory, "data_reduced.Rds")
tse_path <- file.path(data_directory, "tse.Rds")
stage1_linear_03_path <- file.path(
  analysis_directory, "stage1_gender", "tables", "linear_models_0.3M.csv"
)
stage1_linear_36_path <- file.path(
  analysis_directory, "stage1_gender", "tables", "linear_models_3.6M.csv"
)
stage2_wilcoxon_path <- file.path(
  analysis_directory, "stage2", "tables", "wilcoxon_hormones.csv"
)
stage2_growth_path <- file.path(
  analysis_directory, "stage2", "tables",
  "linear_growth_vs_hormones_0.3M.csv"
)
species_alpha_path <- file.path(
  analysis_directory, "stage3", "species", "tables", "alpha_categorical.csv"
)
b_longum_results_path <- file.path(
  analysis_directory, "differential_abundance", "species_clr",
  "early_life_ab", "all_results.tsv"
)
b_longum_values_path <- file.path(
  analysis_directory, "differential_abundance", "species_clr",
  "early_life_ab", "features", "data_transformed.tsv"
)
b_longum_labour_results_path <- file.path(
  analysis_directory, "differential_abundance", "species",
  "Labour.antibiotic", "all_results.tsv"
)
b_longum_labour_values_path <- file.path(
  analysis_directory, "differential_abundance", "species",
  "Labour.antibiotic", "features", "data_transformed.tsv"
)
varibaculum_results_path <- file.path(
  analysis_directory, "differential_abundance", "genus",
  "P.LH_log", "all_results.tsv"
)
varibaculum_values_path <- file.path(
  analysis_directory, "differential_abundance", "genus",
  "P.LH_log", "features", "data_transformed.tsv"
)

required_files <- c(
  gender_data_path, reduced_data_path, tse_path, stage1_linear_03_path,
  stage1_linear_36_path, stage2_wilcoxon_path, stage2_growth_path,
  species_alpha_path, b_longum_results_path,
  b_longum_values_path, b_longum_labour_results_path,
  b_longum_labour_values_path, varibaculum_results_path,
  varibaculum_values_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Required data or analysis output missing: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

# Load the three analysis datasets.
gender_data <- readRDS(gender_data_path)
reduced_data <- readRDS(reduced_data_path)
tse <- readRDS(tse_path)
species_data <- SingleCellExperiment::altExp(tse, "species")
species_metadata <- as.data.frame(SummarizedExperiment::colData(species_data))

# Publication theme and colours.
colour_no <- "#4D4D4D"
colour_yes <- "#2F6B9A"
colour_female <- "#6B6B6B"
colour_male <- "#2F6B9A"

theme_article <- function(base_size = 9) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      legend.position = "none",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

format_p <- function(value) {
  if (is.na(value)) return("NA")
  if (value < 0.001) return(format(value, scientific = TRUE, digits = 2))
  formatC(value, format = "f", digits = 3)
}

save_figure <- function(plot, filename, width, height) {
  stem <- file.path(output_directory, filename)
  ggplot2::ggsave(
    paste0(stem, ".pdf"), plot = plot, width = width, height = height,
    units = "in", device = grDevices::pdf
  )
  ggplot2::ggsave(
    paste0(stem, ".png"), plot = plot, width = width, height = height,
    units = "in", dpi = 600, bg = "white"
  )
  ggplot2::ggsave(
    paste0(stem, ".tiff"), plot = plot, width = width, height = height,
    units = "in", dpi = 600, compression = "lzw", bg = "white"
  )
  invisible(stem)
}

read_feature_values <- function(path, feature, value_name) {
  table <- readr::read_tsv(path, show_col_types = FALSE)
  if (!feature %in% names(table)) {
    stop("Feature not found in MaAsLin3 transformed data: ", feature)
  }
  data.frame(
    Subject_ID = as.character(table[[1]]),
    value = as.numeric(table[[feature]]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::rename(!!value_name := value)
}

# Figure 1A: growth during the first three months by sex.
figure1a_data <- gender_data |>
  dplyr::select(
    Gender, Growth.velocity.0.3M, Weight.gain.0.3M
  ) |>
  tidyr::pivot_longer(
    cols = -Gender, names_to = "Outcome", values_to = "Value"
  ) |>
  dplyr::mutate(
    Gender = factor(Gender, levels = c("female", "male"),
                    labels = c("Female", "Male")),
    Outcome = factor(
      Outcome,
      levels = c("Growth.velocity.0.3M", "Weight.gain.0.3M"),
      labels = c("Growth velocity", "Weight gain")
    )
  )

figure1a <- ggplot2::ggplot(
  figure1a_data, ggplot2::aes(x = Gender, y = Value, fill = Gender)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  ggplot2::geom_jitter(width = 0.13, alpha = 0.28, size = 0.75) +
  ggplot2::facet_wrap(~Outcome, scales = "free_y") +
  ggplot2::scale_fill_manual(values = c(Female = colour_female, Male = colour_male)) +
  ggplot2::labs(x = NULL, y = NULL) +
  theme_article()

# Figure 1B: standardised associations across both growth periods.
stage1_linear <- dplyr::bind_rows(
  readr::read_csv(stage1_linear_03_path, show_col_types = FALSE),
  readr::read_csv(stage1_linear_36_path, show_col_types = FALSE)
)
perinatal_predictors <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy",
  "Birth.height", "Birth.weight"
)
predictor_labels <- c(
  Mothers.age = "Maternal age",
  Pregnancy.weeks = "Pregnancy duration",
  Mother.BMI.before.pregnancy = "Pre-pregnancy BMI",
  BMI.after.pregnancy = "Post-pregnancy BMI",
  Weight.gain.during.pregnancy = "Gestational weight gain",
  Birth.height = "Birth length",
  Birth.weight = "Birth weight"
)
outcome_labels <- c(
  Growth.velocity.0.3M = "Growth velocity\n0-3 months",
  Weight.gain.0.3M = "Weight gain\n0-3 months",
  Growth.velocity.3.6M = "Growth velocity\n3-6 months",
  Weight.gain.3.6M = "Weight gain\n3-6 months"
)

figure1b_data <- stage1_linear |>
  dplyr::filter(predictor %in% perinatal_predictors) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    scale_factor = stats::sd(gender_data[[predictor]], na.rm = TRUE) /
      stats::sd(gender_data[[outcome]], na.rm = TRUE),
    beta = estimate * scale_factor,
    beta_low = (estimate - 1.96 * std.error) * scale_factor,
    beta_high = (estimate + 1.96 * std.error) * scale_factor,
    Predictor = unname(predictor_labels[predictor]),
    Outcome = unname(outcome_labels[outcome]),
    FDR_significant = p_adj_BH < 0.05
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Predictor = factor(Predictor, levels = rev(unname(predictor_labels))),
    Outcome = factor(Outcome, levels = unname(outcome_labels))
  )

figure1b <- ggplot2::ggplot(
  figure1b_data,
  ggplot2::aes(x = beta, y = Predictor, colour = FDR_significant)
) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = beta_low, xmax = beta_high),
    orientation = "y", width = 0.18
  ) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~Outcome, nrow = 1) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = colour_yes)) +
  ggplot2::labs(
    x = "Standardised regression coefficient (95% CI)", y = NULL
  ) +
  theme_article()

# Figure 1C: maternal age and early weight gain.
maternal_age_result <- stage1_linear |>
  dplyr::filter(
    outcome == "Weight.gain.0.3M", predictor == "Mothers.age"
  ) |>
  dplyr::slice(1)

figure1c <- ggplot2::ggplot(
  gender_data,
  ggplot2::aes(x = Mothers.age, y = Weight.gain.0.3M)
) +
  ggplot2::geom_point(alpha = 0.45, size = 1.2, colour = "grey35") +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = colour_yes, fill = "#BFD3E2", linewidth = 0.8
  ) +
  ggplot2::labs(
    x = "Maternal age (years)", y = "Weight gain, birth to 3 months (g)"
  ) +
  theme_article()

# Figure 1D: maternal age and early growth velocity.
maternal_age_growth_result <- stage1_linear |>
  dplyr::filter(
    outcome == "Growth.velocity.0.3M", predictor == "Mothers.age"
  ) |>
  dplyr::slice(1)

figure1d <- ggplot2::ggplot(
  gender_data,
  ggplot2::aes(x = Mothers.age, y = Growth.velocity.0.3M)
) +
  ggplot2::geom_point(alpha = 0.45, size = 1.2, colour = "grey35") +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = colour_yes, fill = "#BFD3E2", linewidth = 0.8
  ) +
  ggplot2::labs(
    x = "Maternal age (years)",
    y = "Growth velocity, birth to 3 months (cm/year)"
  ) +
  theme_article()

figure1 <- patchwork::wrap_plots(
  a = figure1a, b = figure1b, c = figure1c, d = figure1d,
  design = "aa\nbb\ncd"
) +
  patchwork::plot_annotation(tag_levels = "a")

# Figure 2B: pregnancy antibiotic exposure and serum testosterone.
stage2_wilcoxon <- readr::read_csv(stage2_wilcoxon_path, show_col_types = FALSE)
testosterone_result <- stage2_wilcoxon |>
  dplyr::filter(
    response_var == "S.TestoMS_log", var_name == "Pregnancy.antibiotic"
  ) |>
  dplyr::slice(1)

figure2_data <- reduced_data |>
  dplyr::mutate(
    Pregnancy_antibiotic = factor(
      as.character(Pregnancy.antibiotic), levels = c("1", "2"),
      labels = c("No", "Yes")
    )
  )

figure2b <- ggplot2::ggplot(
  figure2_data,
  ggplot2::aes(x = Pregnancy_antibiotic, y = S.TestoMS_log,
               fill = Pregnancy_antibiotic)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.13, alpha = 0.4, size = 1) +
  ggplot2::scale_fill_manual(values = c(No = colour_no, Yes = colour_yes)) +
  ggplot2::labs(
    x = "Maternal antibiotic use during pregnancy",
    y = "Natural log serum testosterone (nmol/L)"
  ) +
  theme_article(base_size = 10)

# Figure 2A: plasma SHBG and weight gain during the first three months.
stage2_growth <- readr::read_csv(stage2_growth_path, show_col_types = FALSE)
shbg_result <- stage2_growth |>
  dplyr::filter(
    outcome == "Weight.gain.0.3M", predictor == "P.SHBG_log"
  ) |>
  dplyr::slice(1)

figure2a <- ggplot2::ggplot(
  reduced_data,
  ggplot2::aes(x = P.SHBG_log, y = Weight.gain.0.3M)
) +
  ggplot2::geom_point(alpha = 0.45, size = 1.2, colour = "grey35") +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = colour_yes, fill = "#BFD3E2", linewidth = 0.8
  ) +
  ggplot2::labs(
    x = "Natural log plasma SHBG (nmol/L)",
    y = "Weight gain, birth to 3 months (g)"
  ) +
  theme_article(base_size = 10)

# Order panels to follow the manuscript: SHBG first, testosterone second.
figure2 <- figure2a + figure2b +
  patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "a")

# Figure 3C: species-level Shannon diversity and early-life antibiotics.
species_alpha <- readr::read_csv(species_alpha_path, show_col_types = FALSE)
shannon_result <- species_alpha |>
  dplyr::filter(variable == "early_life_ab") |>
  dplyr::slice(1)

figure3c_data <- species_metadata |>
  dplyr::mutate(
    Antibiotic_exposure = factor(
      as.character(early_life_ab), levels = c("1", "2"),
      labels = c("No", "Yes")
    )
  )

figure3c <- ggplot2::ggplot(
  figure3c_data,
  ggplot2::aes(x = Antibiotic_exposure, y = shannon_diversity,
               fill = Antibiotic_exposure)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.13, alpha = 0.4, size = 0.9) +
  ggplot2::scale_fill_manual(values = c(No = colour_no, Yes = colour_yes)) +
  ggplot2::labs(
    x = "Early-life antibiotic exposure", y = "Shannon diversity"
  ) +
  theme_article()

# Figure 3B: B. longum CLR abundance and early-life antibiotics.
b_longum_results <- readr::read_tsv(
  b_longum_results_path, show_col_types = FALSE
)
b_longum_result <- b_longum_results |>
  dplyr::filter(
    feature == "Bifidobacterium longum", model == "abundance"
  ) |>
  dplyr::slice(1)
b_longum_values <- read_feature_values(
  b_longum_values_path, "Bifidobacterium longum", "B_longum_clr"
) |>
  dplyr::mutate(
    Antibiotic_exposure = factor(
      as.character(species_metadata$early_life_ab[
        match(Subject_ID, rownames(species_metadata))
      ]),
      levels = c("1", "2"), labels = c("No", "Yes")
    )
  )

figure3b <- ggplot2::ggplot(
  b_longum_values,
  ggplot2::aes(x = Antibiotic_exposure, y = B_longum_clr,
               fill = Antibiotic_exposure)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.13, alpha = 0.4, size = 0.9) +
  ggplot2::scale_fill_manual(values = c(No = colour_no, Yes = colour_yes)) +
  ggplot2::labs(
    x = "Early-life antibiotic exposure",
    y = expression(italic("B. longum") ~ "CLR abundance")
  ) +
  theme_article()

# Alternative Figure 3B: B. longum conditional abundance and labour antibiotics.
# Keep the original CLR/early-life-exposure panel above as the default.
b_longum_labour_results <- readr::read_tsv(
  b_longum_labour_results_path, show_col_types = FALSE
)
b_longum_labour_result <- b_longum_labour_results |>
  dplyr::filter(
    feature == "Bifidobacterium longum", model == "abundance"
  ) |>
  dplyr::slice(1)
b_longum_labour_values <- read_feature_values(
  b_longum_labour_values_path, "Bifidobacterium longum",
  "B_longum_abundance"
) |>
  dplyr::mutate(
    Labour_antibiotic = factor(
      as.character(species_metadata$Labour.antibiotic[
        match(Subject_ID, rownames(species_metadata))
      ]),
      levels = c("1", "2"), labels = c("No", "Yes")
    )
  ) |>
  dplyr::filter(is.finite(B_longum_abundance), !is.na(Labour_antibiotic))

figure3b_labour_option <- ggplot2::ggplot(
  b_longum_labour_values,
  ggplot2::aes(x = Labour_antibiotic, y = B_longum_abundance,
               fill = Labour_antibiotic)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  ggplot2::geom_jitter(width = 0.13, alpha = 0.4, size = 0.9) +
  ggplot2::scale_fill_manual(values = c(No = colour_no, Yes = colour_yes)) +
  ggplot2::labs(
    x = "Labour antibiotic exposure",
    y = expression(italic("B. longum") ~ "conditional abundance")
  ) +
  theme_article()

# Figure 3A: Varibaculum conditional abundance and plasma LH.
varibaculum_results <- readr::read_tsv(
  varibaculum_results_path, show_col_types = FALSE
)
varibaculum_result <- varibaculum_results |>
  dplyr::filter(feature == "Varibaculum", model == "abundance") |>
  dplyr::slice(1)
varibaculum_values <- read_feature_values(
  varibaculum_values_path, "Varibaculum", "Varibaculum_abundance"
) |>
  dplyr::mutate(
    P_LH_log = species_metadata$P.LH_log[
      match(Subject_ID, rownames(species_metadata))
    ]
  ) |>
  dplyr::filter(is.finite(Varibaculum_abundance), is.finite(P_LH_log))

figure3a <- ggplot2::ggplot(
  varibaculum_values,
  ggplot2::aes(x = P_LH_log, y = Varibaculum_abundance)
) +
  ggplot2::geom_point(alpha = 0.65, size = 1.4, colour = "grey30") +
  ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = colour_yes, fill = "#BFD3E2", linewidth = 0.8
  ) +
  ggplot2::labs(
    x = "Natural log plasma LH (U/L)",
    y = expression(italic("Varibaculum") ~ "log relative abundance")
  ) +
  theme_article()

# Order panels to follow the manuscript: Varibaculum, B. longum, diversity.
figure3 <- figure3a + figure3b + figure3c +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(tag_levels = "a")

figure3_labour_option <- figure3a + figure3b_labour_option + figure3c +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(tag_levels = "a")

# Save panel descriptions and statistics separately from the clean figures.
figure_captions <- data.frame(
  Figure = c("1", "1", "1", "1", "2", "2", "3", "3", "3"),
  Panel = c("a", "b", "c", "d", "b", "a", "c", "b", "a"),
  Description = c(
    "Growth velocity and weight gain from birth to 3 months by infant sex.",
    "Standardised linear-regression coefficients for associations between perinatal factors and growth during the 0–3- and 3–6-month periods; blue indicates BH/FDR-adjusted p < 0.05.",
    paste0(
      "Association between maternal age and weight gain from birth to 3 months (estimate ",
      formatC(maternal_age_result$estimate, digits = 3, format = "f"),
      " g per maternal year; BH p = ",
      format_p(maternal_age_result$p_adj_BH), ")."
    ),
    paste0(
      "Association between maternal age and growth velocity from birth to 3 months (estimate ",
      formatC(maternal_age_growth_result$estimate, digits = 3, format = "f"),
      " cm/year per maternal year; BH p = ",
      format_p(maternal_age_growth_result$p_adj_BH), ")."
    ),
    paste0(
      "Serum testosterone by maternal antibiotic use during pregnancy (Mann-Whitney raw p = ",
      format_p(testosterone_result$raw_p), "; BH p = ",
      format_p(testosterone_result$p_adj_BH), ")."
    ),
    paste0(
      "Association between natural-log-transformed plasma SHBG and weight gain from birth to 3 months (estimate ",
      formatC(shbg_result$estimate, digits = 3, format = "f"),
      " g per log-unit; BH p = ", format_p(shbg_result$p_adj_BH), ")."
    ),
    paste0(
      "Species-level Shannon diversity by early-life antibiotic exposure (raw p = ",
      format_p(shannon_result$p.value), "; BH p = ",
      format_p(shannon_result$p.adj), ")."
    ),
    paste0(
      "Bifidobacterium longum CLR abundance by early-life antibiotic exposure (MaAsLin3 raw p = ",
      format_p(b_longum_result$pval_joint), "; q = ",
      format_p(b_longum_result$qval_joint), ")."
    ),
    paste0(
      "Plasma LH and conditional Varibaculum abundance among non-zero samples (n = ",
      nrow(varibaculum_values), "; MaAsLin3 q = ",
      format_p(varibaculum_result$qval_joint), ")."
    )
  ),
  stringsAsFactors = FALSE
)
figure_captions <- figure_captions[order(figure_captions$Figure, figure_captions$Panel), ]
readr::write_csv(
  figure_captions,
  file.path(output_directory, "Figure_captions.csv")
)
readr::write_csv(
  data.frame(
    Figure = "3 (labour-antibiotic option)",
    Panel = "b",
    Description = paste0(
      "Bifidobacterium longum conditional abundance by labour antibiotic ",
      "exposure (MaAsLin3 pval_individual = ",
      format_p(b_longum_labour_result$pval_individual),
      "; qval_individual = ",
      format_p(b_longum_labour_result$qval_individual),
      "; pval_joint = ", format_p(b_longum_labour_result$pval_joint),
      "; qval_joint = ", format_p(b_longum_labour_result$qval_joint),
      "; nominal association)."
    ),
    stringsAsFactors = FALSE
  ),
  file.path(output_directory, "Figure_3B_labour_antibiotic_option_caption.csv")
)

# Save individual panels and complete figures.
save_figure(figure1a, "Figure_1A", 6.5, 3.3)
save_figure(figure1b, "Figure_1B", 11.0, 4.6)
save_figure(figure1c, "Figure_1C", 5.5, 3.3)
save_figure(figure1d, "Figure_1D", 5.5, 3.3)
save_figure(figure1, "Figure_1", 11.0, 10.5)

save_figure(figure2b, "Figure_2B", 4.2, 4.0)
save_figure(figure2a, "Figure_2A", 4.2, 4.0)
save_figure(figure2, "Figure_2", 8.4, 4.0)

# Submission-order copies use the updated panel letters.
save_figure(figure2a, "Figure_2_01_panel_A", 4.2, 4.0)
save_figure(figure2b, "Figure_2_02_panel_B", 4.2, 4.0)

save_figure(figure3c, "Figure_3C", 3.4, 3.4)
save_figure(figure3b, "Figure_3B", 3.4, 3.4)
save_figure(
  figure3b_labour_option, "Figure_3B_labour_antibiotic_option", 3.4, 3.4
)
save_figure(figure3a, "Figure_3A", 3.4, 3.4)
save_figure(figure3, "Figure_3", 10.0, 3.6)
save_figure(figure3a, "Figure_3_01_panel_A", 3.4, 3.4)
save_figure(figure3b, "Figure_3_02_panel_B", 3.4, 3.4)
save_figure(figure3c, "Figure_3_03_panel_C", 3.4, 3.4)
save_figure(
  figure3b_labour_option,
  "Figure_3_02_panel_B_labour_antibiotic_option", 3.4, 3.4
)
save_figure(
  figure3_labour_option, "Figure_3_labour_antibiotic_option", 10.0, 3.6
)

message(
  "Article figures created: ",
  normalizePath(output_directory, mustWork = TRUE)
)
