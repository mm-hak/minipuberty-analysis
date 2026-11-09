# Shared helpers for Stage 3 taxonomic analyses.

# Extract PERMANOVA and PERMDISP results from one RDA model.
extract_rda_table <- function(experiment, name, variable, distance_label) {
  significance <- attr(SingleCellExperiment::reducedDim(experiment, name), "significance")
  extract_one <- function(component, test_type) {
    if (is.null(component)) return(data.frame())
    tab <- as.data.frame(component)
    tab$Term <- rownames(tab)
    tab$Variable <- variable
    tab$Model <- name
    tab$Distance <- distance_label
    tab$TestType <- test_type
    p_column <- intersect(names(tab), c("Pr(>F)", "Pr..F.", "p.value", "p.value..perm.", "P"))
    tab$p_value <- if (length(p_column)) suppressWarnings(as.numeric(tab[[p_column[[1]]]])) else NA_real_
    tab
  }
  permanova <- tryCatch(extract_one(significance$permanova, "PERMANOVA"), error = function(e) data.frame())
  homogeneity <- tryCatch({
    if (is.null(significance$homogeneity)) data.frame()
    else extract_one(significance$homogeneity$anova, "PERMDISP")
  }, error = function(e) data.frame())
  dplyr::bind_rows(permanova, homogeneity) |>
    dplyr::filter(!Term %in% c("Model", "Residual", "Total"))
}

# Run the same RDA analysis for a list of variables or interactions.
run_rda_grid <- function(experiment, variables, assay_type, distance, label,
                         permutations = 999L, prefix = "univariate") {
  results <- list()
  for (variable in variables) {
    name <- safe_file_name(paste("RDA", prefix, label, variable, sep = "_"))
    formula <- stats::as.formula(paste("x ~", variable))
    required_variables <- setdiff(all.vars(formula), "x")
    if (!all(required_variables %in% colnames(SummarizedExperiment::colData(experiment)))) next
    fitted <- tryCatch({
      mia::addRDA(experiment, assay.type = assay_type, formula = formula,
                  distance = distance, na.action = stats::na.exclude,
                  permutations = permutations, name = name)
    }, error = function(e) NULL)
    if (is.null(fitted)) {
      warning("RDA failed for ", variable, " (", label, ")", call. = FALSE)
      next
    }
    experiment <- fitted
    results[[variable]] <- extract_rda_table(experiment, name, variable, label)
  }
  table <- dplyr::bind_rows(results, .id = "id")
  table$p_FDR <- NA_real_
  permanova_rows <- which(table$TestType == "PERMANOVA" & !is.na(table$p_value))
  table$p_FDR[permanova_rows] <- stats::p.adjust(table$p_value[permanova_rows], method = "BH")
  list(experiment = experiment, results = table)
}

# Compare Shannon diversity between categorical groups.
alpha_categorical <- function(experiment, variables, index = "shannon_diversity") {
  metadata <- as.data.frame(SummarizedExperiment::colData(experiment))
  rows <- list()
  for (variable in variables) {
    groups <- sort(unique(stats::na.omit(as.character(metadata[[variable]]))))
    if (length(groups) < 2) next
    for (pair in utils::combn(groups, 2, simplify = FALSE)) {
      first <- metadata[[index]][as.character(metadata[[variable]]) == pair[[1]]]
      second <- metadata[[index]][as.character(metadata[[variable]]) == pair[[2]]]
      first <- stats::na.omit(first); second <- stats::na.omit(second)
      test <- stats::wilcox.test(first, second, exact = FALSE)
      n1 <- length(first); n2 <- length(second); w <- unname(test$statistic)
      rows[[length(rows) + 1L]] <- data.frame(
        variable, group1 = pair[[1]], group2 = pair[[2]], n1, n2,
        mean1 = mean(first), mean2 = mean(second), median1 = median(first),
        median2 = median(second), r_rb = 1 - 2 * (w - n1 * (n1 + 1) / 2) / (n1 * n2),
        p.value = test$p.value
      )
    }
  }
  dplyr::bind_rows(rows) |> dplyr::mutate(p.adj = stats::p.adjust(p.value, "BH"))
}

# Correlate Shannon diversity with continuous variables.
alpha_continuous <- function(experiment, variables, index = "shannon_diversity") {
  metadata <- as.data.frame(SummarizedExperiment::colData(experiment))
  purrr::map_dfr(variables, function(variable) {
    complete <- stats::complete.cases(metadata[, c(index, variable)])
    test <- stats::cor.test(metadata[[index]][complete], metadata[[variable]][complete], method = "spearman")
    data.frame(variable, estimate = unname(test$estimate), p.value = test$p.value, n = sum(complete))
  }) |> dplyr::mutate(p.adj = stats::p.adjust(p.value, "BH"))
}

# Read existing MaAsLin result tables without refitting the models.
read_maaslin_results <- function(maaslin_root, taxa_folder, variables, alpha = 0.05,
                                 interaction = FALSE, interesting_features = character()) {
  purrr::map_dfr(variables, function(variable) {
    interaction_folder <- gsub("[^[:alnum:]_]", "_", variable)
    folder <- if (interaction) c(taxa_folder, "interaction", interaction_folder) else c(taxa_folder, variable)
    path <- do.call(file.path, c(list(maaslin_root), as.list(folder), "all_results.tsv"))
    if (!file.exists(path)) {
      warning("MaAsLin result not found: ", path, call. = FALSE)
      return(data.frame())
    }
    result <- readr::read_tsv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    ) |>
      dplyr::mutate(
        # MaAsLin output contains mixed identifier and error-message types.
        # Read them as text, then restore only statistical columns to numeric.
        dplyr::across(dplyr::any_of(c(
          "coef", "stderr", "pval_individual", "qval_individual",
          "pval_joint", "qval_joint"
        )), ~suppressWarnings(as.numeric(.x)))
      )
    significant <- result |>
      dplyr::filter((!is.na(qval_joint) & qval_joint < alpha) |
                    (!is.na(pval_joint) & pval_joint < alpha))
    keep <- union(significant$feature, interesting_features)
    result |> dplyr::filter(feature %in% keep) |> dplyr::mutate(source_variable = variable)
  })
}

# Plot the most abundant taxa at the selected taxonomic level.
plot_top_composition <- function(experiment, level, top_n = 5L) {
  aggregated <- SingleCellExperiment::altExp(experiment, level)
  abundances <- SummarizedExperiment::assay(aggregated, "relabundance")
  keep <- names(sort(rowSums(abundances, na.rm = TRUE), decreasing = TRUE))[seq_len(min(top_n, nrow(abundances)))]
  as.data.frame(t(abundances[keep, , drop = FALSE])) |>
    tibble::rownames_to_column("sample") |>
    tidyr::pivot_longer(-sample, names_to = "taxon", values_to = "relative_abundance") |>
    ggplot2::ggplot(ggplot2::aes(sample, relative_abundance, fill = taxon)) +
    ggplot2::geom_col() + ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 5)) +
    ggplot2::labs(title = paste("Top", level, "taxa"), x = NULL, y = "Relative abundance")
}
