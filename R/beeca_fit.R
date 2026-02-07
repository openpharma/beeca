#' Quick beeca analysis pipeline
#'
#' @importFrom stats binomial
#'
#' @description
#'
#' A convenience function that streamlines the workflow for conducting a
#' covariate-adjusted marginal treatment effect analysis. This function combines
#' model fitting and marginal effect estimation in a single call, with
#' automatic data preprocessing and informative messages.
#'
#' @details
#'
#' This function provides a simplified interface to the beeca workflow by:
#' - Automatically converting the treatment variable to a factor if needed
#' - Building the model formula from variable names
#' - Fitting the logistic regression model
#' - Computing marginal effects with robust variance estimation
#' - Providing informative progress messages
#'
#' For more control over the analysis, use \code{glm()} followed by
#' \code{get_marginal_effect()} directly.
#'
#' @param data a data frame containing the analysis variables
#' @param outcome character string specifying the outcome variable name.
#' Must be coded as 0/1 or a factor with two levels.
#' @param treatment character string specifying the treatment variable name.
#' Will be converted to a factor if not already one.
#' @param covariates optional character vector specifying covariate names for
#' adjustment. If NULL, an unadjusted analysis is performed.
#' @param method variance estimation method. One of "Ye" (default) or "Ge".
#' See \link[beeca]{get_marginal_effect} for details.
#' @param contrast type of summary measure. One of "diff" (risk difference,
#' default), "rr" (risk ratio), "or" (odds ratio), "logrr" (log risk ratio),
#' or "logor" (log odds ratio).
#' @param reference optional character string or vector specifying reference
#' treatment level(s) for comparisons. If NULL, defaults to first level.
#' @param strata optional character vector specifying stratification variables
#' (only used with method = "Ye").
#' @param family a GLM family. Default is binomial() for logistic regression.
#' @param verbose logical indicating whether to print progress messages.
#' Default is TRUE.
#' @param ... additional arguments passed to \link[beeca]{get_marginal_effect}
#'
#' @return a beeca object (augmented glm object) with marginal effect estimates.
#' See \link[beeca]{get_marginal_effect} for details on the returned object.
#'
#' @export
#'
#' @examples
#' # Simple two-arm analysis
#' fit <- beeca_fit(
#'   data = trial01,
#'   outcome = "aval",
#'   treatment = "trtp",
#'   covariates = "bl_cov",
#'   method = "Ye",
#'   contrast = "diff",
#'   verbose = FALSE
#' )
#'
#' # View results
#' print(fit)
#' summary(fit)
#'
#' @seealso [get_marginal_effect()] for the underlying estimation function
#' @seealso [predict_counterfactuals()] for counterfactual prediction step
#' @seealso [average_predictions()] for averaging step
#' @seealso [estimate_varcov()] for variance estimation step
#' @seealso [apply_contrast()] for contrast computation step
#' @seealso [tidy.beeca()] for extracting results
#' @seealso [summary.beeca()] for detailed output
#' @seealso [print.beeca()] for concise output
#' @seealso [plot.beeca()] and [plot_forest()] for visualizations
#' @seealso [augment.beeca()] for augmented predictions
#' @seealso [as_gt()] for publication-ready tables
beeca_fit <- function(data,
                      outcome,
                      treatment,
                      covariates = NULL,
                      method = c("Ye", "Ge"),
                      contrast = c("diff", "rr", "or", "logrr", "logor"),
                      reference = NULL,
                      strata = NULL,
                      family = binomial(),
                      verbose = TRUE,
                      ...) {

  # Match arguments
  method <- match.arg(method)
  contrast <- match.arg(contrast)

  # Validate inputs
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame", call. = FALSE)
  }

  if (!outcome %in% names(data)) {
    stop(sprintf("Outcome variable '%s' not found in data.\nAvailable variables: %s",
                 outcome,
                 paste(names(data), collapse = ", ")),
         call. = FALSE)
  }

  if (!treatment %in% names(data)) {
    stop(sprintf("Treatment variable '%s' not found in data.\nAvailable variables: %s",
                 treatment,
                 paste(names(data), collapse = ", ")),
         call. = FALSE)
  }

  if (!is.null(covariates)) {
    missing_covs <- covariates[!covariates %in% names(data)]
    if (length(missing_covs) > 0) {
      stop(sprintf("Covariate(s) not found in data: %s",
                   paste(missing_covs, collapse = ", ")),
           call. = FALSE)
    }
  }

  if (!is.null(strata)) {
    missing_strata <- strata[!strata %in% names(data)]
    if (length(missing_strata) > 0) {
      stop(sprintf("Stratification variable(s) not found in data: %s",
                   paste(missing_strata, collapse = ", ")),
           call. = FALSE)
    }
  }

  # Convert treatment to factor if needed
  if (!is.factor(data[[treatment]])) {
    if (verbose) {
      message(sprintf("Converting '%s' to factor", treatment))
    }
    data[[treatment]] <- as.factor(data[[treatment]])
  }

  # Check for missing data
  analysis_vars <- c(outcome, treatment, covariates)
  complete_cases <- stats::complete.cases(data[, analysis_vars, drop = FALSE])
  n_missing <- sum(!complete_cases)

  if (n_missing > 0) {
    if (verbose) {
      message(sprintf("Warning: %d observation(s) with missing data will be excluded",
                      n_missing))
    }
  }

  # Build formula
  if (is.null(covariates)) {
    formula_str <- sprintf("%s ~ %s", outcome, treatment)
    if (verbose) {
      message("Fitting unadjusted logistic regression model")
    }
  } else {
    cov_str <- paste(covariates, collapse = " + ")
    formula_str <- sprintf("%s ~ %s + %s", outcome, treatment, cov_str)
    if (verbose) {
      message(sprintf("Fitting logistic regression with %d covariate(s)",
                      length(covariates)))
    }
  }

  formula_obj <- stats::as.formula(formula_str)

  # Fit model
  fit <- tryCatch(
    stats::glm(formula_obj, family = family, data = data),
    error = function(e) {
      stop(sprintf("Model fitting failed: %s", e$message), call. = FALSE)
    }
  )

  # Check convergence
  if (!fit$converged) {
    warning("Model did not converge. Results may be unreliable.",
            call. = FALSE)
  }

  # Compute marginal effects
  if (verbose) {
    message(sprintf("Computing marginal effects (method: %s, contrast: %s)",
                    method, contrast))
  }

  # Build get_marginal_effect arguments
  gme_args <- list(
    object = fit,
    trt = treatment,
    strata = strata,
    method = method,
    contrast = contrast
  )

  # Only add reference if not NULL (allows get_marginal_effect to use its default)
  if (!is.null(reference)) {
    gme_args$reference <- reference
  }

  # Add any additional arguments
  gme_args <- c(gme_args, list(...))

  result <- tryCatch(
    do.call(get_marginal_effect, gme_args),
    error = function(e) {
      stop(sprintf("Marginal effect estimation failed: %s", e$message),
           call. = FALSE)
    }
  )

  if (verbose) {
    message("Analysis complete!")
  }

  result
}
