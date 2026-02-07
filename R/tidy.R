#' Tidy method for beeca objects
#'
#' @description
#'
#' Extracts and tidies the marginal treatment effect estimates from a beeca object.
#' This function provides a broom-compatible interface for extracting specific
#' statistics from the Analysis Results Data (ARD) structure.
#'
#' @details
#'
#' The `tidy.beeca()` method extracts key inferential statistics from the
#' `marginal_results` component of a beeca object. By default, it returns
#' the contrast estimates (treatment effects) with their standard errors,
#' test statistics, and p-values. Optionally, it can also include the
#' marginal risk estimates for each treatment arm.
#'
#' The function computes Wald test statistics (estimate / std.error) and
#' corresponding two-sided p-values for each estimate.
#'
#' @param x a beeca object (glm object modified by \link[beeca]{get_marginal_effect})
#' @param conf.int logical indicating whether to include confidence intervals.
#' Defaults to FALSE.
#' @param conf.level confidence level for intervals. Defaults to 0.95.
#' @param include_marginal logical indicating whether to include marginal risk
#' estimates for each treatment arm in addition to contrasts. Defaults to FALSE.
#' @param ... additional arguments (not currently used)
#'
#' @return a tibble with columns:
#' \tabular{ll}{
#'  term      \tab The parameter name (contrast or treatment level) \cr
#'  estimate  \tab The point estimate \cr
#'  std.error \tab The standard error of the estimate \cr
#'  statistic \tab The Wald test statistic (estimate / std.error) \cr
#'  p.value   \tab Two-sided p-value from the Wald test \cr
#'  conf.low  \tab Lower confidence limit (if conf.int = TRUE) \cr
#'  conf.high \tab Upper confidence limit (if conf.int = TRUE) \cr
#' }
#'
#' @export
#'
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#' # Tidy the contrast results
#' tidy(fit1)
#'
#' # Include confidence intervals
#' tidy(fit1, conf.int = TRUE)
#'
#' # Include marginal risk estimates for each arm
#' tidy(fit1, include_marginal = TRUE, conf.int = TRUE)
#'
#' @importFrom generics tidy
#' @export
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [print.beeca()] for concise output
#' @seealso [summary.beeca()] for detailed summary output
#' @seealso [augment.beeca()] for augmented predictions
#' @seealso [plot.beeca()] and [plot_forest()] for visualizations
#' @seealso [as_gt()] for publication-ready tables
tidy.beeca <- function(x, conf.int = FALSE, conf.level = 0.95,
                       include_marginal = FALSE, ...) {

  if (!inherits(x, "beeca")) {
    stop("x must be a beeca object (from get_marginal_effect())")
  }

  if (!is.logical(conf.int)) {
    stop("conf.int must be logical (TRUE or FALSE)")
  }

  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a number between 0 and 1")
  }

  if (!is.logical(include_marginal)) {
    stop("include_marginal must be logical (TRUE or FALSE)")
  }

  # Extract marginal_results ARD
  ard <- x$marginal_results

  # Filter for inferential statistics only
  inferential <- ard[ard$ANALTYP1 == "INFERENTIAL", ]

  # Identify contrast statistics (not "risk" or "risk_se")
  contrast_stats <- inferential[!inferential$STAT %in% c("risk", "risk_se"), ]

  # Reshape contrast results: pair estimates with standard errors
  # Get unique contrasts (those not ending in _se)
  contrast_names <- unique(contrast_stats$STAT[!grepl("_se$", contrast_stats$STAT)])

  contrast_results <- lapply(contrast_names, function(contrast_name) {
    se_name <- paste0(contrast_name, "_se")

    # Get rows for this contrast
    est_row <- contrast_stats[contrast_stats$STAT == contrast_name, ]
    se_row <- contrast_stats[contrast_stats$STAT == se_name, ]

    # Match by TRTVAL to handle multiple comparisons
    matched <- merge(
      est_row[, c("TRTVAL", "STATVAL")],
      se_row[, c("TRTVAL", "STATVAL")],
      by = "TRTVAL",
      suffixes = c("_est", "_se")
    )

    data.frame(
      term = matched$TRTVAL,
      estimate = matched$STATVAL_est,
      std.error = matched$STATVAL_se,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, contrast_results)

  # Optionally include marginal risk estimates
  if (include_marginal) {
    marginal_risk <- inferential[inferential$STAT %in% c("risk", "risk_se"), ]

    # Get unique treatment levels
    trt_levels <- unique(marginal_risk$TRTVAL[marginal_risk$STAT == "risk"])

    marginal_results <- lapply(trt_levels, function(trt_level) {
      est_row <- marginal_risk[marginal_risk$TRTVAL == trt_level &
                                 marginal_risk$STAT == "risk", ]
      se_row <- marginal_risk[marginal_risk$TRTVAL == trt_level &
                                marginal_risk$STAT == "risk_se", ]

      data.frame(
        term = paste0("risk_", trt_level),
        estimate = est_row$STATVAL,
        std.error = se_row$STATVAL,
        stringsAsFactors = FALSE
      )
    })

    marginal_df <- do.call(rbind, marginal_results)
    result <- rbind(marginal_df, result)
  }

  # Calculate test statistics and p-values
  result$statistic <- result$estimate / result$std.error
  result$p.value <- 2 * stats::pnorm(-abs(result$statistic))

  # Add confidence intervals if requested
  if (conf.int) {
    alpha <- 1 - conf.level
    z_crit <- stats::qnorm(1 - alpha / 2)

    result$conf.low <- result$estimate - z_crit * result$std.error
    result$conf.high <- result$estimate + z_crit * result$std.error
  }

  # Convert to tibble and return
  dplyr::as_tibble(result)
}
