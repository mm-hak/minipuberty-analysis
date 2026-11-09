# Shared helpers for the publication analysis scripts.

check_and_load_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Required package(s) are not installed: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(packages, function(package) {
    suppressPackageStartupMessages(
      library(package, character.only = TRUE, quietly = TRUE)
    )
  }))
}

ensure_output_directories <- function(paths) {
  invisible(lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE))
}

safe_file_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

save_table <- function(x, path) {
  readr::write_csv(as.data.frame(x), path, na = "")
  invisible(path)
}

save_plot <- function(plot, stem, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = paste0(stem, ".pdf"), plot = plot,
    width = width, height = height, units = "in"
  )
  ggplot2::ggsave(
    filename = paste0(stem, ".png"), plot = plot,
    width = width, height = height, units = "in", dpi = 300
  )
  invisible(c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png")))
}

plot_histograms_by_group <- function(df, outcomes, group_var) {
  df |>
    dplyr::select(dplyr::all_of(c(group_var, outcomes))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(outcomes),
      names_to = "Variable", values_to = "Value"
    ) |>
    dplyr::filter(!is.na(.data[[group_var]]), !is.na(Value)) |>
    ggplot2::ggplot(ggplot2::aes(x = Value, fill = as.factor(.data[[group_var]]))) +
    ggplot2::geom_histogram(position = "dodge", alpha = 0.7, bins = 30) +
    ggplot2::facet_wrap(~Variable, scales = "free") +
    ggplot2::scale_fill_manual(values = c("darkblue", "darkgray")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Value", y = "Count", fill = group_var)
}

plot_histograms_by_group_with_gender <- function(
    df, outcomes, group_var, gender_var = "Gender") {
  df |>
    dplyr::select(dplyr::all_of(c(group_var, gender_var, outcomes))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(outcomes),
      names_to = "Variable", values_to = "Value"
    ) |>
    dplyr::filter(
      !is.na(.data[[group_var]]), !is.na(.data[[gender_var]]), !is.na(Value)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = Value, fill = as.factor(.data[[group_var]]))) +
    ggplot2::geom_histogram(position = "dodge", alpha = 0.7, bins = 30) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data[[gender_var]]),
      cols = ggplot2::vars(Variable), scales = "free"
    ) +
    ggplot2::scale_fill_manual(values = c("darkblue", "darkgray")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Value", y = "Count", fill = group_var)
}

plot_pair <- function(predictor, response, data, gender = FALSE) {
  mapping <- if (gender) {
    ggplot2::aes(
      x = .data[[predictor]], y = .data[[response]],
      color = Gender, shape = Gender
    )
  } else {
    ggplot2::aes(x = .data[[predictor]], y = .data[[response]])
  }

  ggplot2::ggplot(data, mapping) +
    ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 2) +
    ggplot2::geom_smooth(
      method = "lm", formula = y ~ x, se = FALSE, linetype = "dashed"
    ) +
    ggplot2::labs(x = predictor, y = response) +
    ggplot2::theme_minimal()
}

make_class_distribution <- function(data, variables, by_gender = FALSE) {
  long <- data |>
    dplyr::select(dplyr::all_of(c(variables, if (by_gender) "Gender"))) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(variables), names_to = "variable", values_to = "level"
    ) |>
    dplyr::mutate(level = as.character(level))

  grouping <- c(if (by_gender) "Gender", "variable", "level")
  long |>
    dplyr::count(dplyr::across(dplyr::all_of(grouping)), name = "n") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(if (by_gender) "Gender", "variable")))) |>
    dplyr::mutate(percent = 100 * n / sum(n)) |>
    dplyr::ungroup()
}

make_wilcox_table <- function(data, outcomes, predictors) {
  results <- list()
  failures <- list()
  index <- 1L
  failure_index <- 1L

  for (outcome in outcomes) {
    for (predictor in predictors) {
      subset <- data |>
        dplyr::select(dplyr::all_of(c(outcome, predictor))) |>
        stats::na.omit()

      test_result <- tryCatch({
        test <- rstatix::wilcox_test(
          subset, formula = stats::as.formula(paste0(outcome, " ~ ", predictor))
        )
        effect <- rstatix::wilcox_effsize(
          subset, formula = stats::as.formula(paste0(outcome, " ~ ", predictor))
        )
        levels <- sort(unique(as.character(subset[[predictor]])))
        if (length(levels) != 2) stop("Predictor does not have exactly two levels.")

        value1 <- subset[[outcome]][subset[[predictor]] == levels[[1]]]
        value2 <- subset[[outcome]][subset[[predictor]] == levels[[2]]]
        data.frame(
          response_var = outcome, var_name = predictor,
          group1 = levels[[1]], group2 = levels[[2]],
          mean1 = mean(value1, na.rm = TRUE), mean2 = mean(value2, na.rm = TRUE),
          median1 = stats::median(value1, na.rm = TRUE),
          median2 = stats::median(value2, na.rm = TRUE),
          U = test$statistic[[1]], raw_p = test$p[[1]],
          effect_size_r = effect$effsize[[1]]
        )
      }, error = function(error) {
        failures[[failure_index]] <<- data.frame(
          response_var = outcome,
          var_name = predictor,
          error = conditionMessage(error)
        )
        failure_index <<- failure_index + 1L
        NULL
      })

      if (!is.null(test_result)) {
        results[[index]] <- test_result
        index <- index + 1L
      }
    }
  }

  result <- dplyr::bind_rows(results)
  if (nrow(result) == 0) {
    failure_text <- paste(
      dplyr::bind_rows(failures)$error,
      collapse = "; "
    )
    stop("All Wilcoxon analyses failed: ", failure_text, call. = FALSE)
  }
  if (length(failures) > 0) {
    warning(
      length(failures),
      " Wilcoxon analysis/analyses failed; inspect the input levels.",
      call. = FALSE
    )
  }
  result$p_adj_BH <- stats::p.adjust(result$raw_p, method = "BH")
  result
}

make_significant_boxplots <- function(results, data, top_n = 3) {
  selected <- results[results$p_adj_BH < 0.05, , drop = FALSE]
  subtitle <- "Significant differences (FDR < 0.05)"
  if (nrow(selected) == 0) {
    selected <- results[order(-abs(results$effect_size_r)), , drop = FALSE]
    selected <- utils::head(selected, top_n)
    subtitle <- "Top effect-size differences (not FDR-significant)"
  }

  plots <- lapply(seq_len(nrow(selected)), function(i) {
    row <- selected[i, ]
    ggplot2::ggplot(
      data, ggplot2::aes(x = .data[[row$var_name]], y = .data[[row$response_var]])
    ) +
      ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[row$var_name]])) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(
        x = row$var_name, y = row$response_var,
        title = paste(row$response_var, "by", row$var_name), subtitle = subtitle
      )
  })
  patchwork::wrap_plots(plots, ncol = 1)
}

run_art_models <- function(data, outcomes, predictors, time_period) {
  results <- list()
  failures <- list()
  index <- 1L
  failure_index <- 1L

  for (predictor in predictors) {
    for (outcome in outcomes) {
      formula <- stats::as.formula(sprintf("%s ~ %s * Gender", outcome, predictor))
      result <- tryCatch({
        model <- ARTool::art(formula, data = data)
        table <- as.data.frame(stats::anova(model))
        residual_df <- table$Df.res[[1]]
        table |>
          tibble::rownames_to_column("Effect") |>
          dplyr::transmute(
            Time = time_period, Outcome = outcome, Predictor = predictor,
            Effect, F_stat = `F value`, raw_p = `Pr(>F)`,
            effect_size = `F value` / (`F value` + residual_df)
          )
      }, error = function(error) {
        failures[[failure_index]] <<- data.frame(
          Time = time_period, Outcome = outcome, Predictor = predictor,
          error = conditionMessage(error)
        )
        failure_index <<- failure_index + 1L
        NULL
      })
      if (!is.null(result)) {
        results[[index]] <- result
        index <- index + 1L
      }
    }
  }

  result_table <- dplyr::bind_rows(results) |>
    dplyr::mutate(
      p_adj_BH = stats::p.adjust(raw_p, method = "BH"),
      effect_magnitude = dplyr::case_when(
        effect_size < 0.01 ~ "negligible",
        effect_size < 0.06 ~ "small",
        effect_size < 0.14 ~ "medium",
        TRUE ~ "large"
      )
    )
  list(results = result_table, failures = dplyr::bind_rows(failures))
}

fit_linear_model_grid <- function(data, outcomes, predictors, gender_interaction = FALSE) {
  models <- list()
  tables <- list()
  index <- 1L

  for (outcome in outcomes) {
    for (predictor in predictors) {
      formula_text <- if (gender_interaction) {
        sprintf("%s ~ %s * Gender", outcome, predictor)
      } else {
        sprintf("%s ~ %s", outcome, predictor)
      }
      model_id <- safe_file_name(paste(outcome, predictor, sep = "__"))
      model <- stats::lm(stats::as.formula(formula_text), data = data)
      models[[model_id]] <- model
      summary <- summary(model)
      tables[[index]] <- broom::tidy(model) |>
        dplyr::filter(term != "(Intercept)") |>
        dplyr::transmute(
          model = model_id, outcome = outcome, predictor = predictor,
          term, estimate, std.error, statistic,
          r_squared = summary$r.squared, p.value
        )
      index <- index + 1L
    }
  }

  results <- dplyr::bind_rows(tables) |>
    dplyr::mutate(p_adj_BH = stats::p.adjust(p.value, method = "BH"))
  list(results = results, models = models)
}

save_model_diagnostics <- function(models, directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- character(length(models))
  names(paths) <- names(models)
  for (name in names(models)) {
    path <- file.path(directory, paste0(name, ".png"))
    grDevices::png(path, width = 1200, height = 1200, res = 150)
    old_parameters <- graphics::par(mfrow = c(2, 2))
    graphics::plot(models[[name]])
    graphics::par(old_parameters)
    grDevices::dev.off()
    paths[[name]] <- path
  }
  paths
}

plot_model_estimates <- function(results, top_n = 10, interaction_only = FALSE) {
  selected <- results
  if (interaction_only) selected <- selected[grepl(":", selected$term), , drop = FALSE]
  significant <- selected[selected$p_adj_BH < 0.05, , drop = FALSE]
  if (nrow(significant) > 0) {
    selected <- significant
    subtitle <- "FDR-significant estimates"
  } else {
    selected <- selected[order(selected$p.value), , drop = FALSE]
    selected <- utils::head(selected, top_n)
    subtitle <- "Lowest unadjusted p-values (none FDR-significant)"
  }
  selected <- selected |>
    dplyr::mutate(
      label = paste(outcome, term, sep = " | "),
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )

  ggplot2::ggplot(selected, ggplot2::aes(x = estimate, y = reorder(label, estimate))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = conf.low, xmax = conf.high),
      orientation = "y", width = 0.2
    ) +
    ggplot2::geom_point(color = "darkblue") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Estimate (95% CI)", y = NULL, subtitle = subtitle)
}
