#' Convert beeca marginal_results to cards ARD format
#'
#' @description
#' Converts the `marginal_results` component from a beeca analysis to the
#' ARD (Analysis Results Data) format used by the cards package. This enables
#' integration with the cards/cardx ecosystem for comprehensive reporting and
#' quality control workflows.
#'
#' @param marginal_results A beeca `marginal_results` tibble/data frame
#'   containing the output from `get_marginal_effect()$marginal_results`.
#'
#' @return A cards ARD object (tibble with class 'card') containing:
#'   * `group1`: Treatment variable name (character)
#'   * `group1_level`: Treatment level (list-column, numeric if possible, else character)
#'   * `variable`: Outcome variable name (character)
#'   * `variable_level`: Specific level of variable (character, NA for beeca)
#'   * `stat_name`, `stat_label`: Statistic identifier and human-readable label
#'   * `stat`: The calculated value (list-column)
#'   * `context`: Analysis context (combining ANALTYP1 and ANALMETH)
#'   * `fmt_fn`, `warning`, `error`: Cards-specific metadata (list-columns)
#'
#' @details
#' The function maps beeca's CDISC-inspired ARD structure to the cards package
#' format:
#'
#' | beeca | cards | Notes |
#' |-------|-------|-------|
#' | `TRTVAR` | `group1` | Treatment variable name |
#' | `TRTVAL` | `group1_level` | Treatment level |
#' | `PARAM` | `variable` | Outcome variable |
#' | `STAT` | `stat_name` | Statistic identifier |
#' | `STATVAL` | `stat` | Value (numeric vs list-column) |
#' | `ANALTYP1` + `ANALMETH` | `context` | Analysis context |
#'
#' Original beeca metadata is preserved in the `beeca_description` attribute.
#'
#' @section Package Requirements:
#' This function requires the `cards` package to be installed. It is listed
#' as a suggested dependency and will provide an informative error if not
#' available.
#'
#' @examples
#' if (requireNamespace("cards", quietly = TRUE)) {
#'   # Fit model and get beeca results
#'   trial01$trtp <- factor(trial01$trtp)
#'   fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'     get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#'   # Convert to cards format
#'   cards_ard <- beeca_to_cards_ard(fit1$marginal_results)
#'
#'   # Print the cards ARD (uses print method for 'card' class)
#'   print(cards_ard)
#'
#'   # Bind with other cards ARDs
#'   combined_ard <- cards::bind_ard(
#'     cards_ard,
#'     cards::ard_continuous(trial01, by = trtp, variables = bl_cov)
#'   )
#' }
#'
#' @seealso [get_marginal_effect()] for creating the input `marginal_results`
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [as_gt()] for publication-ready tables
#' @seealso [tidy.beeca()] for broom-compatible output
#' @seealso `vignette("ard-cards-integration")` for detailed integration examples
#'
#' @export
beeca_to_cards_ard <- function(marginal_results) {

  if (!requireNamespace("cards", quietly = TRUE)) {
    stop("Package 'cards' is required for this conversion. Please install it with:\n  install.packages(\"cards\")")
  }

  # Create stat_label mapping
  stat_label_map <- c(
    "N" = "N",
    "n" = "n",
    "%" = "%",
    "risk" = "Risk",
    "risk_se" = "Risk SE",
    "diff" = "Difference",
    "diff_se" = "Difference SE",
    "rr" = "Risk Ratio",
    "rr_se" = "Risk Ratio SE",
    "or" = "Odds Ratio",
    "or_se" = "Odds Ratio SE",
    "logrr" = "Log Risk Ratio",
    "logrr_se" = "Log Risk Ratio SE",
    "logor" = "Log Odds Ratio",
    "logor_se" = "Log Odds Ratio SE"
  )

  result <- marginal_results |>
    dplyr::rename(
      group1 = .data$TRTVAR,
      variable = .data$PARAM,
      stat_name = .data$STAT
    ) |>
    dplyr::mutate(
      # Convert group1_level to list-column (required by cards for compatibility)
      # Try to convert to numeric if possible, otherwise keep as character
      group1_level = lapply(.data$TRTVAL, function(x) {
        num_val <- suppressWarnings(as.numeric(x))
        if (!is.na(num_val)) num_val else x
      }),

      # Add human-readable labels
      stat_label = dplyr::case_when(
        .data$stat_name %in% names(stat_label_map) ~ stat_label_map[.data$stat_name],
        TRUE ~ toupper(.data$stat_name)
      ),

      # Consolidate context from multiple beeca columns
      context = paste(
        tolower(.data$ANALTYP1),
        .data$ANALMETH,
        sep = "_"
      ),

      # Convert atomic numeric to list column (required by cards)
      stat = as.list(.data$STATVAL),

      # Add empty list columns for cards compatibility
      fmt_fn = list(NULL),
      warning = list(NULL),
      error = list(NULL),

      # Add variable_level (not used in beeca's current structure)
      variable_level = NA_character_
    ) |>
    dplyr::select(
      .data$group1, .data$group1_level, .data$variable, .data$variable_level,
      .data$stat_name, .data$stat_label, .data$stat, .data$context,
      .data$fmt_fn, .data$warning, .data$error
    )

  # Store beeca metadata as attributes
  attr(result, "beeca_description") <- unique(marginal_results[["ANALDESC"]])

  # Apply cards column ordering and class
  result |>
    cards::tidy_ard_column_order() |>
    cards::as_card()
}
