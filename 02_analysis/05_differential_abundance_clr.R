# MaAsLin3 differential-abundance analysis using CLR-transformed abundances.
# Run via workflow.R or with 02_analysis as the working directory.

# Load helper functions and packages ------------------------------------------
source("analysis_helpers.R", local = FALSE)
source("differential_abundance_helpers.R", local = FALSE)
load_differential_abundance_packages()

# Load the processed microbiome data ------------------------------------------
input_path <- file.path("..", "data", "tse.Rds")
if (!file.exists(input_path)) stop("Processed input not found: ", input_path, call. = FALSE)
tse <- readRDS(input_path)

# Define variables excluded from CLR models ----------------------------------
# Use the original growth and weight-gain variables, matching the relative-
# abundance, alpha-diversity, and Bray-Curtis analyses. Log-transformed
# hormone variables remain because those are the prespecified hormone scales.
drop_variables <- c(
  "dbp_dominance", "faith_diversity", "observed_richness", "shannon_diversity",
  "Growth.velocity.0.3M_log", "Growth.velocity.3.6M_log",
  "Weight.gain.0.3M_log", "Weight.gain.3.6M_log",
  "Growth.velocity.0.6M", "Growth.velocity.0.6M_log",
  "Weight.gain.0.6M", "Weight.gain.0.6M_log",
  "S.TestoMS", "P.LH", "P.FSH", "P.SHBG"
)

# Define interaction models ---------------------------------------------------
interactions <- c(
  "Labour.antibiotic * Intervention",
  "Pregnancy.antibiotic * Intervention",
  "Labour.antibiotic * Weight.gain.during.pregnancy",
  "Pregnancy.antibiotic * Weight.gain.during.pregnancy",
  "Labour.antibiotic * Mother.BMI.before.pregnancy",
  "Pregnancy.antibiotic * Mother.BMI.before.pregnancy",
  "Intervention * Weight.gain.during.pregnancy"
)

# Run species, genus, and phylum models ---------------------------------------
output_root <- file.path("output", "differential_abundance")
for (taxa_level in c("species", "genus", "phylum")) {
  message("Running CLR MaAsLin3 models: ", taxa_level)
  run_maaslin_level(
    tse = tse,
    taxa_level = taxa_level,
    output_root = output_root,
    drop_variables = drop_variables,
    interactions = interactions,
    assay = "clr"
  )
}

message("CLR differential analysis complete: ", output_root)
