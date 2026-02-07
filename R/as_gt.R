#' Convert beeca object to gt table
#'
#' @description
#'
#' Creates a publication-ready clinical trial table from beeca analysis results
#' using the gt package. The table includes marginal risks by treatment arm,
#' treatment effect estimates with confidence intervals, and supports
#' customization for titles, footnotes, and analysis set information.
#'
#' @param x a beeca object (from \link[beeca]{get_marginal_effect} or \link[beeca]{beeca_fit})
#' @param title character string for table title. Default is NULL (no title).
#' @param subtitle character string for table subtitle. Default is NULL.
#' @param source_note character string for table source/footnote. Default is NULL.
#' @param analysis_set character string describing the analysis population
#'   (e.g., "Full Analysis Set (FAS)", "Per Protocol Set"). Default is NULL.
#' @param analysis_set_n integer, number of subjects in analysis set.
#'   If NULL (default), uses the sample size from the model.
#' @param conf.level numeric, confidence level for intervals. Default is 0.95.
#' @param risk_digits integer, decimal places for risk estimates. Default is 1
#'   (displays as percentages).
#' @param effect_digits integer, decimal places for treatment effect. Default is 2.
#' @param include_ci logical, whether to include confidence intervals. Default is TRUE.
#' @param include_pvalue logical, whether to include p-values. Default is TRUE.
#' @param risk_percent logical, if TRUE displays risks as percentages (0-100),
#'   if FALSE displays as proportions (0-1). Default is TRUE.
#' @param ... additional arguments passed to gt functions
#'
#' @return a gt table object
#'
#' @export
#'
#' @examples
#' if (requireNamespace("gt", quietly = TRUE)) {
#'   # Fit model
#'   trial01$trtp <- factor(trial01$trtp)
#'   fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'     get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#'   # Create clinical trial table
#'   as_gt(fit,
#'     title = "Table 14.2.1: Primary Efficacy Analysis",
#'     subtitle = "Response Rate by Treatment Group",
#'     source_note = paste("Risk difference estimated using g-computation",
#'                         "with robust variance (Ye et al. 2023)"),
#'     analysis_set = "Full Analysis Set (FAS)"
#'   )
#' }
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [tidy.beeca()] for data frame output
#' @seealso [summary.beeca()] for console output
#' @seealso [print.beeca()] for concise output
#' @seealso [plot.beeca()] and [plot_forest()] for visualizations
#' @seealso [beeca_summary_table()] for simpler data frame alternative
#' @seealso [beeca_to_cards_ard()] for cards ARD integration
as_gt <- function(x, ...) {
  UseMethod("as_gt")
}

#' @export
as_gt.default <- function(x, ...) {
  stop("x must be a beeca object (from get_marginal_effect() or beeca_fit())")
}

#' @rdname as_gt
#' @export
as_gt.beeca <- function(x,
                        title = NULL,
                        subtitle = NULL,
                        source_note = NULL,
                        analysis_set = NULL,
                        analysis_set_n = NULL,
                        conf.level = 0.95,
                        risk_digits = 1,
                        effect_digits = 2,
                        include_ci = TRUE,
                        include_pvalue = TRUE,
                        risk_percent = TRUE,
                        ...) {

  if (!requireNamespace("gt", quietly = TRUE)) {
    stop("Package 'gt' is required for this function. Please install it with install.packages('gt').")
  }

  # Validation
 if (!inherits(x, "beeca")) {
    stop("x must be a beeca object (from get_marginal_effect() or beeca_fit())")
  }

  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a number between 0 and 1")
  }

  # Extract information from beeca object
  ard <- x$marginal_results
  n_obs <- if (is.null(analysis_set_n)) nrow(stats::model.frame(x)) else analysis_set_n

  # Get variance method and contrast type for footnotes
  variance_method <- attr(x$marginal_se, "type")
  if (is.null(variance_method)) variance_method <- "robust"

  contrast_attr <- attr(x$marginal_est, "contrast")
  contrast_type <- if (!is.null(contrast_attr)) {
    sub(":.*", "", contrast_attr[1])
  } else {
    "diff"
  }

  # Get contrast label
  contrast_labels <- c(
    "diff" = "Risk Difference",
    "rr" = "Risk Ratio",
    "or" = "Odds Ratio",
    "logrr" = "Log Risk Ratio",
    "logor" = "Log Odds Ratio"
  )
  contrast_label <- contrast_labels[contrast_type]
  if (is.na(contrast_label)) contrast_label <- contrast_type

  # Get reference arm
  reference_arm <- attr(x$marginal_est, "reference")
  if (is.null(reference_arm)) reference_arm <- names(x$counterfactual.means)[1]

  # Build marginal risks table
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  ci_label <- sprintf("%d%% CI", round(conf.level * 100))

  # Extract per-arm statistics from ARD
  trt_levels <- unique(ard$TRTVAL[ard$ANALTYP1 == "DESCRIPTIVE"])

  arm_data <- lapply(trt_levels, function(trt) {
    desc <- ard[ard$TRTVAL == trt & ard$ANALTYP1 == "DESCRIPTIVE", ]
    inf <- ard[ard$TRTVAL == trt & ard$ANALTYP1 == "INFERENTIAL", ]

    n_total <- desc$STATVAL[desc$STAT == "N"]
    n_resp <- desc$STATVAL[desc$STAT == "n"]
    pct <- desc$STATVAL[desc$STAT == "%"]
    risk <- inf$STATVAL[inf$STAT == "risk"]
    risk_se <- inf$STATVAL[inf$STAT == "risk_se"]

    data.frame(
      treatment = trt,
      n = n_total,
      responders = n_resp,
      observed_pct = pct,
      risk = risk,
      risk_se = risk_se,
      stringsAsFactors = FALSE
    )
  })

  arm_df <- do.call(rbind, arm_data)

  # Format risk estimates
  multiplier <- if (risk_percent) 100 else 1
  arm_df$risk_fmt <- sprintf(paste0("%.", risk_digits, "f"), arm_df$risk * multiplier)
  arm_df$risk_ci_low <- arm_df$risk - z_crit * arm_df$risk_se
  arm_df$risk_ci_high <- arm_df$risk + z_crit * arm_df$risk_se

  arm_df$risk_ci_fmt <- sprintf(
    paste0("(%.", risk_digits, "f, %.", risk_digits, "f)"),
    arm_df$risk_ci_low * multiplier,
    arm_df$risk_ci_high * multiplier
  )

  # Build treatment effects data
  effects_data <- lapply(seq_along(x$marginal_est), function(i) {
    est <- x$marginal_est[i]
    se <- x$marginal_se[i]
    comparison <- attr(x$marginal_est, "contrast")[i]

    ci_low <- est - z_crit * se
    ci_high <- est + z_crit * se
    z_stat <- est / se
    p_val <- 2 * stats::pnorm(-abs(z_stat))

    # Format based on contrast type
    if (contrast_type %in% c("diff")) {
      est_fmt <- sprintf(paste0("%.", effect_digits, "f"), est * multiplier)
      ci_fmt <- sprintf(
        paste0("(%.", effect_digits, "f, %.", effect_digits, "f)"),
        ci_low * multiplier, ci_high * multiplier
      )
    } else {
      est_fmt <- sprintf(paste0("%.", effect_digits, "f"), est)
      ci_fmt <- sprintf(
        paste0("(%.", effect_digits, "f, %.", effect_digits, "f)"),
        ci_low, ci_high
      )
    }

    data.frame(
      comparison = comparison,
      estimate = est,
      estimate_fmt = est_fmt,
      se = se,
      ci_low = ci_low,
      ci_high = ci_high,
      ci_fmt = ci_fmt,
      z_stat = z_stat,
      p_value = p_val,
      p_value_fmt = format_pvalue(p_val),
      stringsAsFactors = FALSE
    )
  })

  effects_df <- do.call(rbind, effects_data)

  # Create main table data
  table_rows <- list()

  # Add header for marginal risks section
  for (i in seq_len(nrow(arm_df))) {
    row <- arm_df[i, ]
    is_ref <- row$treatment == reference_arm

    row_data <- data.frame(
      category = "Marginal Risk",
      label = row$treatment,
      n = as.character(row$n),
      responders = as.character(row$responders),
      estimate = row$risk_fmt,
      ci = row$risk_ci_fmt,
      p_value = if (is_ref) "" else NA_character_,
      stringsAsFactors = FALSE
    )
    table_rows[[length(table_rows) + 1]] <- row_data
  }

  # Add treatment effects section
  for (i in seq_len(nrow(effects_df))) {
    row <- effects_df[i, ]
    row_data <- data.frame(
      category = contrast_label,
      label = row$comparison,
      n = "",
      responders = "",
      estimate = row$estimate_fmt,
      ci = row$ci_fmt,
      p_value = row$p_value_fmt,
      stringsAsFactors = FALSE
    )
    table_rows[[length(table_rows) + 1]] <- row_data
  }

  table_df <- do.call(rbind, table_rows)

  # Build gt table
  gt_table <- gt::gt(table_df, groupname_col = "category") |>
    gt::cols_label(
      label = "Treatment",
      n = "N",
      responders = "Responders",
      estimate = if (risk_percent) "Estimate (%)" else "Estimate",
      ci = ci_label,
      p_value = "P-value"
    )

  # Hide columns based on options
  if (!include_ci) {
    gt_table <- gt_table |> gt::cols_hide(columns = "ci")
  }

  if (!include_pvalue) {
    gt_table <- gt_table |> gt::cols_hide(columns = "p_value")
  }

  # Add title and subtitle
  if (!is.null(title) || !is.null(subtitle)) {
    gt_table <- gt_table |>
      gt::tab_header(
        title = title,
        subtitle = subtitle
      )
  }

  # Add analysis set info
  if (!is.null(analysis_set)) {
    analysis_set_text <- if (!is.null(analysis_set_n)) {
      sprintf("%s (N = %d)", analysis_set, analysis_set_n)
    } else {
      sprintf("%s (N = %d)", analysis_set, n_obs)
    }

    gt_table <- gt_table |>
      gt::tab_source_note(source_note = analysis_set_text)
  }

  # Add source note / footnote
  if (!is.null(source_note)) {
    gt_table <- gt_table |>
      gt::tab_source_note(source_note = source_note)
  }

  # Add method footnote
  method_note <- sprintf(
    "Marginal risks estimated using g-computation. Variance estimated using %s method.",
    variance_method
  )
  gt_table <- gt_table |>
    gt::tab_source_note(source_note = method_note)

  # Style the table
  gt_table <- gt_table |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_row_groups()
    ) |>
    gt::tab_style(
      style = gt::cell_text(align = "right"),
      locations = gt::cells_body(columns = c("n", "responders", "estimate", "ci", "p_value"))
    ) |>
    gt::tab_style(
      style = gt::cell_text(align = "right"),
      locations = gt::cells_column_labels(columns = c("n", "responders", "estimate", "ci", "p_value"))
    ) |>
    gt::tab_options(
      table.font.size = gt::px(12),
      heading.title.font.size = gt::px(14),
      heading.subtitle.font.size = gt::px(12),
      source_notes.font.size = gt::px(10)
    )

  gt_table
}


#' Format p-value for clinical trial tables
#'
#' @description
#'
#' Formats p-values according to common clinical trial reporting conventions.
#' Values less than 0.001 are displayed as "<0.001", otherwise rounded to
#' 3 decimal places.
#'
#' @param p numeric p-value(s) to format
#' @param digits integer, number of decimal places. Default is 3.
#' @param small_threshold numeric, values below this threshold are displayed
#'   as "<threshold". Default is 0.001.
#'
#' @return character vector of formatted p-values
#'
#' @export
#'
#' @examples
#' format_pvalue(0.0001)
#' format_pvalue(c(0.05, 0.001, 0.0001))
format_pvalue <- function(p, digits = 3, small_threshold = 0.001) {
  sapply(p, function(pval) {
    if (is.na(pval)) {
      return(NA_character_)
    }
    if (pval < small_threshold) {
      return(sprintf("<%s", format(small_threshold, scientific = FALSE)))
    }
    format(round(pval, digits), nsmall = digits)
  })
}


#' Create a summary statistics table for beeca analysis
#'
#' @description
#'
#' Creates a summary table showing descriptive statistics and treatment
#' effects from a beeca analysis. This function provides a simpler alternative
#' to `as_gt()` that returns a data frame suitable for further customization.
#'
#' @param x a beeca object
#' @param conf.level numeric, confidence level for intervals. Default is 0.95.
#' @param risk_percent logical, if TRUE displays risks as percentages. Default is TRUE.
#'
#' @return a data frame with summary statistics
#'
#' @export
#'
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#' beeca_summary_table(fit)
beeca_summary_table <- function(x, conf.level = 0.95, risk_percent = TRUE) {

  if (!inherits(x, "beeca")) {
    stop("x must be a beeca object")
  }

  ard <- x$marginal_results
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  multiplier <- if (risk_percent) 100 else 1

  # Get treatment levels
  trt_levels <- unique(ard$TRTVAL[ard$ANALTYP1 == "DESCRIPTIVE"])

  # Build per-arm summary
  arm_summary <- lapply(trt_levels, function(trt) {
    desc <- ard[ard$TRTVAL == trt & ard$ANALTYP1 == "DESCRIPTIVE", ]
    inf <- ard[ard$TRTVAL == trt & ard$ANALTYP1 == "INFERENTIAL", ]

    risk <- inf$STATVAL[inf$STAT == "risk"]
    risk_se <- inf$STATVAL[inf$STAT == "risk_se"]

    data.frame(
      treatment = trt,
      N = desc$STATVAL[desc$STAT == "N"],
      n_responders = desc$STATVAL[desc$STAT == "n"],
      observed_rate = desc$STATVAL[desc$STAT == "%"],
      marginal_risk = risk * multiplier,
      marginal_risk_se = risk_se * multiplier,
      marginal_risk_ci_low = (risk - z_crit * risk_se) * multiplier,
      marginal_risk_ci_high = (risk + z_crit * risk_se) * multiplier,
      stringsAsFactors = FALSE
    )
  })

  arm_df <- do.call(rbind, arm_summary)

  # Get contrast info
  contrast_type <- sub(":.*", "", attr(x$marginal_est, "contrast")[1])

  # Build treatment effects summary
  effects_summary <- data.frame(
    comparison = attr(x$marginal_est, "contrast"),
    estimate = x$marginal_est,
    std_error = x$marginal_se,
    stringsAsFactors = FALSE
  )

  effects_summary$ci_low <- effects_summary$estimate - z_crit * effects_summary$std_error
  effects_summary$ci_high <- effects_summary$estimate + z_crit * effects_summary$std_error
  effects_summary$z_statistic <- effects_summary$estimate / effects_summary$std_error
  effects_summary$p_value <- 2 * stats::pnorm(-abs(effects_summary$z_statistic))

  # Scale if risk difference and percent requested
  if (contrast_type == "diff" && risk_percent) {
    effects_summary$estimate <- effects_summary$estimate * 100
    effects_summary$std_error <- effects_summary$std_error * 100
    effects_summary$ci_low <- effects_summary$ci_low * 100
    effects_summary$ci_high <- effects_summary$ci_high * 100
  }

  row.names(effects_summary) <- NULL

  list(
    arm_statistics = dplyr::as_tibble(arm_df),
    treatment_effects = dplyr::as_tibble(effects_summary),
    metadata = list(
      conf_level = conf.level,
      contrast_type = contrast_type,
      variance_method = attr(x$marginal_se, "type"),
      reference = attr(x$marginal_est, "reference"),
      n_total = nrow(stats::model.frame(x)),
      risk_percent = risk_percent
    )
  )
}
