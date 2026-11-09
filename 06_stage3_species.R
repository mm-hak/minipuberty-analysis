# Stage 3 species-level microbiome analysis.
# Run via workflow.R or with 02_analysis as the working directory.

# Load helper functions and packages.
source("analysis_helpers.R", local = FALSE)
source("stage3_helpers.R", local = FALSE)
check_and_load_packages(c(
  "tidyverse", "broom", "mia", "miaViz", "patchwork",
  "SummarizedExperiment", "SingleCellExperiment"
))

# Load processed microbiome data.
taxa_level <- "species"
input_path <- file.path("..", "data", "tse.Rds")
if (!file.exists(input_path)) {
  stop("Processed input not found: ", input_path, call. = FALSE)
}
tse <- readRDS(input_path)
data <- SingleCellExperiment::altExp(tse, taxa_level)

# Define output folders.
output_root <- file.path("output", "stage3", taxa_level)
paths <- list(
  tables = file.path(output_root, "tables"),
  figures = file.path(output_root, "figures"),
  objects = file.path(output_root, "objects")
)
ensure_output_directories(paths)

# Define analysis variables.
categorical <- c(
  "Pregnancy.antibiotic", "Labour.antibiotic",
  "Child.antibiotic.under.7.days", "Child.antibiotic.7.days.to.3M",
  "early_life_ab", "Intervention", "Birth.mode"
)
continuous <- c(
  "Mothers.age", "Pregnancy.weeks", "Mother.BMI.before.pregnancy",
  "BMI.after.pregnancy", "Weight.gain.during.pregnancy",
  "Birth.height", "Birth.weight", "Growth.velocity.0.3M",
  "Growth.velocity.3.6M", "Weight.gain.0.3M", "Weight.gain.3.6M",
  "P.LH_log", "P.FSH_log", "P.SHBG_log"
)
outcomes <- c(
  "Growth.velocity.0.3M", "Growth.velocity.3.6M",
  "Weight.gain.0.3M", "Weight.gain.3.6M"
)
drop <- c(
  "dbp_dominance", "faith_diversity", "observed_richness",
  "shannon_diversity", "P.LH", "P.FSH", "P.SHBG", "S.TestoMS",
  "S.TestoMS_sqrt", "S.Testo_group", paste0(outcomes, "_log"),
  "Growth.velocity.0.6M", "Growth.velocity.0.6M_log",
  "Weight.gain.0.6M", "Weight.gain.0.6M_log"
)
variables <- setdiff(colnames(SummarizedExperiment::colData(tse)), drop)
variables <- variables[make.names(variables) == variables]

# Define the variables used to read CLR-based MaAsLin3 results.
clr_drop <- c(
  "dbp_dominance", "faith_diversity", "observed_richness",
  "shannon_diversity", paste0(outcomes, "_log"),
  "Growth.velocity.0.6M", "Growth.velocity.0.6M_log",
  "Weight.gain.0.6M", "Weight.gain.0.6M_log",
  "S.TestoMS", "P.LH", "P.FSH", "P.SHBG"
)
clr_variables <- setdiff(
  colnames(SummarizedExperiment::colData(tse)), clr_drop
)
clr_variables <- clr_variables[make.names(clr_variables) == clr_variables]

interactions <- c(
  "Labour.antibiotic * Intervention",
  "Pregnancy.antibiotic * Intervention",
  "Labour.antibiotic * Weight.gain.during.pregnancy",
  "Pregnancy.antibiotic * Weight.gain.during.pregnancy",
  "Labour.antibiotic * Mother.BMI.before.pregnancy",
  "Pregnancy.antibiotic * Mother.BMI.before.pregnancy",
  "early_life_ab * Intervention"
)
clr_interactions <- c(
  "Labour.antibiotic * Intervention",
  "Pregnancy.antibiotic * Intervention",
  "Labour.antibiotic * Weight.gain.during.pregnancy",
  "Pregnancy.antibiotic * Weight.gain.during.pregnancy",
  "Labour.antibiotic * Mother.BMI.before.pregnancy",
  "Pregnancy.antibiotic * Mother.BMI.before.pregnancy",
  "Intervention * Weight.gain.during.pregnancy"
)

# Community composition.
composition <- plot_top_composition(tse, taxa_level)
save_plot(
  composition,
  file.path(paths$figures, "community_composition_top_taxa"),
  12, 7
)

# Shannon diversity distribution.
shannon_plot <- ggplot2::ggplot(
  as.data.frame(SummarizedExperiment::colData(data)),
  ggplot2::aes(x = shannon_diversity)
) +
  ggplot2::geom_histogram(
    binwidth = 0.1, fill = "skyblue", color = "black"
  ) +
  ggplot2::geom_density(
    ggplot2::aes(y = after_stat(count) * 0.1),
    alpha = 0.2, fill = "blue"
  ) +
  ggplot2::theme_minimal()
save_plot(
  shannon_plot,
  file.path(paths$figures, "shannon_distribution"),
  8, 6
)

# Alpha diversity analysis.
alpha_cat <- alpha_categorical(data, categorical)
alpha_cont <- alpha_continuous(data, continuous)
save_table(alpha_cat, file.path(paths$tables, "alpha_categorical.csv"))
save_table(alpha_cont, file.path(paths$tables, "alpha_continuous.csv"))

# Bray-Curtis beta diversity analysis.
bray <- run_rda_grid(
  data, variables, "relabundance", "bray", "Bray-Curtis"
)
save_table(
  bray$results,
  file.path(paths$tables, "beta_bray_univariate.csv")
)

# Aitchison beta diversity analysis.
aitchison <- run_rda_grid(
  bray$experiment, variables, "clr", "euclidean", "Aitchison"
)
save_table(
  aitchison$results,
  file.path(paths$tables, "beta_aitchison_univariate.csv")
)

# Alpha diversity interaction analysis.
interaction_alpha <- purrr::map_dfr(interactions, function(interaction) {
  model <- stats::lm(
    stats::as.formula(paste("shannon_diversity ~", interaction)),
    data = as.data.frame(SummarizedExperiment::colData(data))
  )
  broom::tidy(model, conf.int = TRUE) |>
    dplyr::mutate(Model = interaction)
}) |>
  dplyr::filter(term != "(Intercept)") |>
  dplyr::mutate(p_adj_BH = stats::p.adjust(p.value, "BH"))
save_table(
  interaction_alpha,
  file.path(paths$tables, "alpha_interactions.csv")
)

# Bray-Curtis interaction analysis.
bray_interactions <- run_rda_grid(
  aitchison$experiment, interactions, "relabundance", "bray",
  "Bray-Curtis", prefix = "interaction"
)
save_table(
  bray_interactions$results,
  file.path(paths$tables, "beta_bray_interactions.csv")
)

# Aitchison interaction analysis.
aitchison_interactions <- run_rda_grid(
  bray_interactions$experiment, interactions, "clr", "euclidean",
  "Aitchison", prefix = "interaction"
)
save_table(
  aitchison_interactions$results,
  file.path(paths$tables, "beta_aitchison_interactions.csv")
)

# Locate differential abundance results.
maaslin_root <- file.path("output", "differential_abundance")
required_maaslin_directories <- file.path(
  maaslin_root, c(taxa_level, paste0(taxa_level, "_clr"))
)
missing_maaslin_directories <- required_maaslin_directories[
  !dir.exists(required_maaslin_directories)
]
if (length(missing_maaslin_directories) > 0) {
  stop(
    "Required differential-abundance output is missing: ",
    paste(missing_maaslin_directories, collapse = ", "),
    ". Run 04_differential_abundance_relative.R and ",
    "05_differential_abundance_clr.R before Stage 3.",
    call. = FALSE
  )
}

# Select significant and prespecified taxa.
interesting_features <- c(
  "Bifidobacterium_longum", "Bifidobacterium_breve",
  "Bifidobacterium_bifidum", "Bifidobacterium_adolescentis",
  "Bifidobacterium_dentium"
)
maaslin_relative <- read_maaslin_results(
  maaslin_root, taxa_level, variables,
  interesting_features = interesting_features
)
maaslin_clr <- read_maaslin_results(
  maaslin_root, paste0(taxa_level, "_clr"), clr_variables,
  interesting_features = interesting_features
)
maaslin_interactions <- read_maaslin_results(
  maaslin_root, taxa_level, interactions,
  interaction = TRUE, interesting_features = interesting_features
)
maaslin_interactions_clr <- read_maaslin_results(
  maaslin_root, paste0(taxa_level, "_clr"), clr_interactions,
  interaction = TRUE, interesting_features = interesting_features
)

# Save differential abundance tables.
save_table(
  maaslin_relative,
  file.path(paths$tables, "differential_abundance_relative.csv")
)
save_table(
  maaslin_clr,
  file.path(paths$tables, "differential_abundance_clr.csv")
)
save_table(
  maaslin_interactions,
  file.path(paths$tables, "differential_abundance_interactions_relative.csv")
)
save_table(
  maaslin_interactions_clr,
  file.path(paths$tables, "differential_abundance_interactions_clr.csv")
)

# Save complete Stage 3 results.
stage3_results <- list(
  alpha = list(
    categorical = alpha_cat,
    continuous = alpha_cont,
    interactions = interaction_alpha
  ),
  beta = list(
    bray = bray$results,
    aitchison = aitchison$results,
    bray_interactions = bray_interactions$results,
    aitchison_interactions = aitchison_interactions$results
  ),
  differential_abundance = list(
    relative = maaslin_relative,
    clr = maaslin_clr,
    interactions_relative = maaslin_interactions,
    interactions_clr = maaslin_interactions_clr
  )
)
saveRDS(
  stage3_results,
  file.path(paths$objects, paste0("stage3_", taxa_level, "_results.rds"))
)

message("Stage 3 ", taxa_level, " analysis complete: ", output_root)
