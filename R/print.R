#' Print method for beeca objects
#'
#' @description
#'
#' Provides a concise summary of a beeca analysis, showing the treatment effect
#' estimate, standard error, and p-value. For more detailed output, use
#' \code{summary()}.
#'
#' @param x a beeca object (from \link[beeca]{get_marginal_effect})
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
#' # Concise output
#' print(fit1)
#'
#' # More detail
#' summary(fit1)
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [summary.beeca()] for detailed output
#' @seealso [tidy.beeca()] for broom-compatible output
#' @seealso [plot.beeca()] for visualizations
#' @seealso [augment.beeca()] for augmented predictions
#' @seealso [as_gt()] for publication-ready tables
print.beeca <- function(x, digits = 4, ...) {

  cat("beeca: Covariate-Adjusted Marginal Treatment Effect\n")
  cat(strrep("=", 55), "\n\n", sep = "")

  # Get contrast information
  contrast_label <- attr(x$marginal_est, "contrast")
  if (is.null(contrast_label)) {
    contrast_label <- "Treatment Effect"
  }

  # Calculate p-value
  z_stat <- x$marginal_est / x$marginal_se
  p_val <- 2 * stats::pnorm(-abs(z_stat))

  # Number of comparisons
  n_contrasts <- length(x$marginal_est)

  if (n_contrasts == 1) {
    # Single contrast
    cat(sprintf("Contrast:  %s\n", contrast_label))
    cat(sprintf("Estimate:  %s (SE = %s)\n",
                format(round(x$marginal_est, digits), nsmall = digits),
                format(round(x$marginal_se, digits), nsmall = digits)))
    cat(sprintf("Z-value:   %s\n", format(round(z_stat, digits), nsmall = digits)))
    cat(sprintf("P-value:   %s\n", format.pval(p_val, digits = digits)))
  } else {
    # Multiple contrasts
    cat(sprintf("%d treatment comparisons:\n\n", n_contrasts))
    for (i in seq_along(x$marginal_est)) {
      cat(sprintf("  %s:\n", contrast_label[i]))
      cat(sprintf("    Estimate = %s (SE = %s), p = %s\n",
                  format(round(x$marginal_est[i], digits), nsmall = digits),
                  format(round(x$marginal_se[i], digits), nsmall = digits),
                  format.pval(p_val[i], digits = digits)))
    }
  }

  cat("\n")
  cat("Use summary() for detailed results\n")
  cat("Use tidy() for broom-compatible output\n")
  cat("Use plot() for visualizations\n")

  invisible(x)
}
