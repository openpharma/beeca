#' Summary method for beeca objects
#'
#' @description
#'
#' Provides a comprehensive summary of a beeca analysis, including marginal
#' risks for each treatment arm, treatment effect estimates with confidence
#' intervals, and key model information.
#'
#' @param object a beeca object (from \link[beeca]{get_marginal_effect})
#' @param conf.level confidence level for confidence intervals. Default is 0.95.
#' @param digits integer, number of decimal places for rounding. Default is 4.
#' @param ... additional arguments (not currently used)
#'
#' @return Invisibly returns the input object
#'
#' @export
#'
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#' # Detailed summary
#' summary(fit1)
#'
#' # With 90% confidence intervals
#' summary(fit1, conf.level = 0.90)
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [print.beeca()] for concise output
#' @seealso [tidy.beeca()] for broom-compatible output
#' @seealso [plot.beeca()] for visualizations
#' @seealso [augment.beeca()] for augmented predictions
#' @seealso [as_gt()] for publication-ready tables
summary.beeca <- function(object, conf.level = 0.95, digits = 4, ...) {

  cat("\n")
  cat("Covariate-Adjusted Marginal Treatment Effect Analysis\n")
  cat(strrep("=", 60), "\n\n", sep = "")

  # Model information
  cat("Model Information:\n")
  cat(strrep("-", 60), "\n", sep = "")

  variance_method <- attr(object$marginal_se, "type")
  if (is.null(variance_method)) {
    variance_method <- "Unknown"
  }

  cat(sprintf("  Working model:      %s\n",
              "Logistic regression with covariate adjustment"))
  cat(sprintf("  Variance estimator: %s\n", variance_method))

  # Extract contrast type from attributes
  contrast_attr <- attr(object$marginal_est, "contrast")
  if (!is.null(contrast_attr)) {
    contrast_type <- sub(":.*", "", contrast_attr[1])
  } else {
    contrast_type <- "Unknown"
  }
  cat(sprintf("  Contrast type:      %s\n", contrast_type))

  # Get sample size from model data
  n_obs <- nrow(stats::model.frame(object))
  cat(sprintf("  Sample size:        %d\n", n_obs))

  # Get outcome variable name
  outcome_var <- as.character(object$formula[[2]])
  cat(sprintf("  Outcome variable:   %s\n", outcome_var))

  cat("\n")

  # Marginal risks table
  cat("Marginal Risks (g-computation):\n")
  cat(strrep("-", 60), "\n", sep = "")

  risks_df <- data.frame(
    Treatment = names(object$counterfactual.means),
    Risk = round(object$counterfactual.means, digits),
    SE = round(sqrt(diag(object$robust_varcov)), digits),
    stringsAsFactors = FALSE
  )

  # Calculate confidence intervals for risks
  z_crit <- stats::qnorm(1 - (1 - conf.level) / 2)
  risks_df$CI_Lower <- round(risks_df$Risk - z_crit * risks_df$SE, digits)
  risks_df$CI_Upper <- round(risks_df$Risk + z_crit * risks_df$SE, digits)

  # Format CI
  risks_df$CI <- sprintf("(%s, %s)",
                         format(risks_df$CI_Lower, nsmall = digits),
                         format(risks_df$CI_Upper, nsmall = digits))

  # Print table
  risks_print <- risks_df[, c("Treatment", "Risk", "SE", "CI")]
  names(risks_print)[4] <- sprintf("%d%% CI", round(conf.level * 100))

  print(risks_print, row.names = FALSE, right = FALSE)
  cat("\n")

  # Treatment effects table
  cat("Treatment Effect Estimates:\n")
  cat(strrep("-", 60), "\n", sep = "")

  effects_df <- data.frame(
    Comparison = attr(object$marginal_est, "contrast"),
    Estimate = round(object$marginal_est, digits),
    SE = round(object$marginal_se, digits),
    stringsAsFactors = FALSE
  )

  # Calculate test statistics and p-values
  effects_df$Z_value <- round(effects_df$Estimate / effects_df$SE, digits)
  effects_df$P_value <- 2 * stats::pnorm(-abs(effects_df$Z_value))

  # Calculate confidence intervals
  effects_df$CI_Lower <- round(effects_df$Estimate - z_crit * effects_df$SE, digits)
  effects_df$CI_Upper <- round(effects_df$Estimate + z_crit * effects_df$SE, digits)

  # Format outputs
  effects_df$CI <- sprintf("(%s, %s)",
                           format(effects_df$CI_Lower, nsmall = digits),
                           format(effects_df$CI_Upper, nsmall = digits))
  effects_df$P_value <- format.pval(effects_df$P_value, digits = max(1, digits - 2))

  # Print table
  effects_print <- effects_df[, c("Comparison", "Estimate", "SE", "Z_value", "CI", "P_value")]
  names(effects_print)[4] <- "Z value"
  names(effects_print)[5] <- sprintf("%d%% CI", round(conf.level * 100))
  names(effects_print)[6] <- "P value"

  print(effects_print, row.names = FALSE, right = FALSE)
  cat("\n")

  # Footer
  cat(strrep("-", 60), "\n", sep = "")
  cat("Note: Standard errors and tests based on robust variance estimation\n")
  cat("Use tidy() for data frame output or plot() for visualizations\n")

  invisible(object)
}
