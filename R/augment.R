#' Augment method for beeca objects
#'
#' @description
#'
#' Augments the original dataset with predictions and counterfactual outcomes
#' from a beeca analysis. This function provides a broom-compatible interface
#' for accessing individual-level predictions and potential outcomes under
#' different treatment assignments.
#'
#' @details
#'
#' The `augment.beeca()` method returns the original dataset used in the
#' analysis, augmented with additional columns containing:
#' - Fitted values from the working model (`.fitted`)
#' - Counterfactual predictions for each treatment level
#' - Optionally, residuals and other diagnostic information
#'
#' This is particularly useful for:
#' - Examining individual predictions and potential outcomes
#' - Creating plots of treatment effects across covariates
#' - Conducting sensitivity analyses
#' - Understanding how g-computation works at the individual level
#'
#' Each counterfactual column represents the predicted outcome if that
#' subject were assigned to a specific treatment level, holding all other
#' covariates constant.
#'
#' @param x a beeca object (glm object modified by \link[beeca]{get_marginal_effect})
#' @param data optional dataset to augment. If NULL (default), uses the data
#' from the original model fit.
#' @param newdata deprecated. Use `data` instead.
#' @param type.predict type of prediction to include as `.fitted`. Options are
#' "response" (default, predicted probabilities) or "link" (linear predictor scale).
#' @param type.residuals type of residuals to include. Options are "deviance" (default),
#' "pearson", "response", or "working". Set to NULL to exclude residuals.
#' @param ... additional arguments (not currently used)
#'
#' @return a tibble containing the original data augmented with:
#' \tabular{ll}{
#'  .rownames \tab Row names from the original data (if present) \cr
#'  .fitted   \tab Fitted values from the working model \cr
#'  .resid    \tab Residuals (if type.residuals is not NULL) \cr
#'  .counterfactual_\[level\] \tab Predicted outcome if assigned to each treatment level \cr
#' }
#'
#' @importFrom generics augment
#' @export
#'
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#' # Augment with counterfactual predictions
#' augmented <- augment(fit1)
#' head(augmented)
#'
#' # Access counterfactual predictions for treatment level 1
#' augmented$.counterfactual_1
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # Examine predictions by baseline covariate
#'   library(ggplot2)
#'   ggplot(augmented, aes(x = bl_cov, y = .counterfactual_1 - .counterfactual_0)) +
#'     geom_point() +
#'     labs(y = "Individual Treatment Effect")
#' }
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [predict_counterfactuals()] for the underlying prediction method
#' @seealso [tidy.beeca()] for tidied parameter estimates
#' @seealso [print.beeca()] for concise output
#' @seealso [summary.beeca()] for detailed summary output
#' @seealso [plot.beeca()] and [plot_forest()] for visualizations
#' @seealso [as_gt()] for publication-ready tables
augment.beeca <- function(x, data = NULL, newdata = NULL,
                          type.predict = "response",
                          type.residuals = "deviance", ...) {

  if (!inherits(x, "beeca")) {
    stop("x must be a beeca object (from get_marginal_effect())")
  }

  # Handle deprecated newdata argument
  if (!is.null(newdata)) {
    warning("'newdata' is deprecated; please use 'data' instead.",
            call. = FALSE)
    if (is.null(data)) {
      data <- newdata
    }
  }

  # Get data from model if not provided
  if (is.null(data)) {
    data <- stats::model.frame(x)
  }

  # Validate type.predict
  if (!type.predict %in% c("response", "link")) {
    stop("type.predict must be 'response' or 'link'")
  }

  # Validate type.residuals
  if (!is.null(type.residuals)) {
    valid_resid_types <- c("deviance", "pearson", "response", "working")
    if (!type.residuals %in% valid_resid_types) {
      stop("type.residuals must be one of: ", paste(valid_resid_types, collapse = ", "),
           " or NULL")
    }
  }

  # Start with the data
  result <- data

  # Add row names if present
  if (!is.null(rownames(data))) {
    result$.rownames <- rownames(data)
  }

  # Add fitted values from the working model
  result$.fitted <- stats::predict(x, type = type.predict, newdata = data)

  # Add residuals if requested
  # Only add residuals if using the original model data
  # (residuals don't make sense for new/different data)
  if (!is.null(type.residuals)) {
    original_data <- stats::model.frame(x)
    if (nrow(data) == nrow(original_data) &&
        isTRUE(all.equal(rownames(data), rownames(original_data)))) {
      result$.resid <- stats::residuals(x, type = type.residuals)
    } else {
      # For custom data, we can't compute residuals from the fit
      # So we skip adding them
      warning("Residuals not added for custom data; only available for original model data",
              call. = FALSE)
    }
  }

  # Add counterfactual predictions for each treatment level
  # Only add counterfactuals if using the original model data
  # (counterfactuals in the beeca object correspond to original data)
  original_data <- stats::model.frame(x)
  if (nrow(data) == nrow(original_data) &&
      isTRUE(all.equal(rownames(data), rownames(original_data)))) {

    # Extract counterfactual predictions from the beeca object
    # This is a tibble where each column is a treatment level
    cf_predictions <- x$counterfactual.predictions

    # Get treatment variable name from attributes
    trt_var <- attr(cf_predictions, "treatment.variable")

    # Get column names (treatment levels)
    trt_levels <- names(cf_predictions)

    # Add each counterfactual prediction as a column
    for (trt_level in trt_levels) {
      col_name <- paste0(".counterfactual_", trt_level)
      result[[col_name]] <- cf_predictions[[trt_level]]
    }
  } else {
    # For custom data, counterfactuals aren't available
    warning("Counterfactual predictions not added for custom data; only available for original model data",
            call. = FALSE)
  }

  # Convert to tibble and return
  dplyr::as_tibble(result)
}
