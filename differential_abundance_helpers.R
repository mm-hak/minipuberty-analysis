# Shared helpers for MaAsLin3 differential-abundance analyses.

# Check and load the packages required for MaAsLin3.
load_differential_abundance_packages <- function() {
  check_and_load_packages(c(
    "maaslin3", "mia", "S4Vectors", "SummarizedExperiment",
    "SingleCellExperiment", "TreeSummarizedExperiment"
  ))
}

# Remove metadata variables that are not included in a model family.
prepare_maaslin_experiment <- function(tse, taxa_level, drop_variables) {
  experiment <- SingleCellExperiment::altExp(tse, taxa_level)
  metadata <- as.data.frame(SummarizedExperiment::colData(experiment))
  metadata <- metadata[, !names(metadata) %in% drop_variables, drop = FALSE]
  SummarizedExperiment::colData(experiment) <- S4Vectors::DataFrame(metadata)

  # Keep taxa above 0.1% relative abundance in at least 10 samples.
  experiment <- mia::subsetByPrevalent(
    experiment,
    assay.type = "relabundance",
    detection = 0.001,
    prevalence = 10 / ncol(experiment)
  )

  # Remove samples that are empty after taxon filtering.
  relative_abundance <- SummarizedExperiment::assay(experiment, "relabundance")
  empty_samples <- which(colSums(relative_abundance) == 0)
  if (length(empty_samples) > 0) {
    message(
      "Removing ", length(empty_samples), " empty sample(s) at ", taxa_level,
      " level: ", paste(colnames(experiment)[empty_samples], collapse = ", ")
    )
    experiment <- experiment[, -empty_samples, drop = FALSE]
  }
  experiment
}

# Run one MaAsLin3 model with the original relative-abundance settings.
run_maaslin_relative_model <- function(experiment, formula, output_directory,
                                       warn_prevalence = FALSE) {
  maaslin3::maaslin3(
    input_data = SummarizedExperiment::assay(experiment, "relabundance"),
    input_metadata = as.data.frame(SummarizedExperiment::colData(experiment)),
    output = output_directory,
    formula = formula,
    normalization = "NONE",
    transform = "LOG",
    zero_threshold = 0,
    warn_prevalence = warn_prevalence,
    standardize = TRUE,
    verbosity = "ERROR",
    plot_associations = TRUE
  )
}

# Run one MaAsLin3 model with the original CLR settings.
run_maaslin_clr_model <- function(experiment, formula, output_directory) {
  maaslin3::maaslin3(
    input_data = SummarizedExperiment::assay(experiment, "clr_filtered"),
    input_metadata = as.data.frame(SummarizedExperiment::colData(experiment)),
    output = output_directory,
    formula = formula,
    normalization = "NONE",
    transform = "NONE",
    zero_threshold = -Inf,
    min_abundance = -Inf,
    warn_prevalence = FALSE,
    standardize = TRUE,
    verbosity = "ERROR",
    plot_associations = TRUE,
    evaluate_only = "abundance"
  )
}

# Run all univariate and interaction models for one taxonomic level.
run_maaslin_level <- function(tse, taxa_level, output_root, drop_variables,
                              interactions, assay = c("relative", "clr")) {
  assay <- match.arg(assay)
  experiment <- prepare_maaslin_experiment(tse, taxa_level, drop_variables)

  if (assay == "clr") {
    experiment <- mia::transformAssay(
      experiment,
      method = "clr",
      assay.type = "relabundance",
      pseudocount = TRUE,
      name = "clr_filtered"
    )
  }

  level_directory <- file.path(
    output_root,
    if (assay == "clr") paste0(taxa_level, "_clr") else taxa_level
  )
  dir.create(level_directory, recursive = TRUE, showWarnings = FALSE)

  # Univariate models ---------------------------------------------------------
  metadata_variables <- colnames(SummarizedExperiment::colData(experiment))
  for (variable in metadata_variables) {
    formula <- stats::as.formula(paste("~", variable))
    model_directory <- file.path(level_directory, variable)
    if (assay == "relative") {
      run_maaslin_relative_model(experiment, formula, model_directory)
    } else {
      run_maaslin_clr_model(experiment, formula, model_directory)
    }
  }

  # Interaction models --------------------------------------------------------
  for (interaction in unique(interactions)) {
    formula <- stats::as.formula(paste("~", interaction))
    safe_interaction <- gsub("[^[:alnum:]_]", "_", interaction)
    model_directory <- file.path(level_directory, "interaction", safe_interaction)
    if (assay == "relative") {
      run_maaslin_relative_model(
        experiment, formula, model_directory, warn_prevalence = TRUE
      )
    } else {
      run_maaslin_clr_model(experiment, formula, model_directory)
    }
  }

  invisible(level_directory)
}
