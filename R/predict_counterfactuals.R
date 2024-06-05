#' Predict counterfactual outcomes in GLM models
#'
#' @description
#'
#' This function calculates counterfactual predictions for each level of a
#' specified treatment variable in a generalized linear model (GLM).
#' It is designed to aid in the assessment of treatment effects by predicting
#' outcomes under different treatments under causal inference framework.
#'
#' @param object a fitted \code{\link[stats]{glm}} object for which
#' counterfactual predictions are desired.
#'
#' @param trt a string specifying the name of the treatment variable
#' in the model formula. It must be one of the linear predictor variables used
#' in fitting the `object`.
#'
#' @return an updated `glm` object appended with an
#' additional component `counterfactual.predictions`.
#'
#' This component contains a tibble with two columns: `cf_pred_0` and `cf_pred_1`,
#' representing counterfactual predictions for each level of the
#' treatment variable. A descriptive `label` attribute explains the counterfactual
#' scenario associated with each column.
#'
#' @importFrom stats predict
#' @importFrom dplyr as_tibble
#' @export
#'
#' @details
#' The function works by creating two new datasets from the original data used
#' to fit the GLM model. In these datasets, the treatment variable
#' is set to each of its levels across all records (e.g., patients).
#'
#' Predictions are then made for each dataset based on the fitted GLM model,
#' simulating the response variable under each treatment condition.
#'
#' The results are stored in a tidy format and appended to the original model
#' object for further analysis or inspection.
#'
#' For averaging counterfactual outcomes, apply `average_predictions()`.
#'
#' @seealso [average_predictions()] for averaging counterfactual
#' predictions.
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#'
#' @examples
#' # Preparing data and fitting a GLM model
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
#'
#' # Generating counterfactual predictions
#' fit2 <- predict_counterfactuals(fit1, "trtp")
#'
#' # Accessing the counterfactual predictions
#' fit2$counterfactual.predictions
#' attributes(fit2$counterfactual.predictions)
#'
predict_counterfactuals <- function(object, trt) {
  object <- .assert_sanitized(object, trt)

  data <- .get_data(object)

  # calculate counterfactual predictions for each trt level
  X0 <- X1 <- data
  X0[, trt] <- levels(data[[trt]])[1]
  X1[, trt] <- levels(data[[trt]])[2]
  cf_pred_0 <- predict(object, newdata = X0, type = "response")
  cf_pred_1 <- predict(object, newdata = X1, type = "response")

  # store predictions in tidy dataframe
  cf_predictions <- data.frame(
    cf_pred_0 = cf_pred_0,
    cf_pred_1 = cf_pred_1
  ) |>
    dplyr::as_tibble()

  cf_lbl <- \(trt, x) paste0("Counterfactual predictions setting ", trt, "=", x, " for all records")
  attr(cf_predictions, "label") <- c(cf_lbl(trt, levels(data[[trt]])[1]), cf_lbl(trt, levels(data[[trt]])[2]))
  attr(cf_predictions, "treatment.variable") <- trt

  object$counterfactual.predictions <- cf_predictions

  return(object)
}
