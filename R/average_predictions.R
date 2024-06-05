#' Average over counterfactual predictions
#'
#' @description
#'
#' `average_predictions()` averages counterfactual predictions stored within
#' a `glm` object. This is pivotal for estimating treatment contrasts and
#' associated variance estimates using g-computation. The function assumes
#' predictions are generated via `predict_counterfactuals()`.
#'
#' @details
#' The `average_predictions()` function calculates the average over the
#' counterfactual predictions which can then be used to estimate a treatment
#' contrast and associated variance estimate.
#'
#' The function appends a `glm` object with the
#' averaged counterfactual predictions.
#'
#' Note: Ensure that the `glm` object has been adequately prepared with
#' `predict_counterfactuals()` before applying `average_predictions()`.
#' Failure to do so may result in errors indicating missing components.
#'
#' @param object a fitted \code{\link[stats]{glm}} object augmented with
#' counterfactual predictions named: `counterfactual.predictions`
#'
#' @return an updated `glm` object appended with an additional component
#' `counterfactual.means`.
#'
#' @export
#'
#' @examples
#'
#' # Use the trial01 dataset
#' data(trial01)
#'
#' # ensure the treatment indicator is a factor
#' trial01$trtp <- factor(trial01$trtp)
#'
#' # fit glm model for trial data
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
#'
#' # Preprocess fit1 as required by average_predictions
#' fit2 <- fit1 |>
#'   predict_counterfactuals(trt = "trtp")
#'
#' # average over the counterfactual predictions
#' fit3 <- average_predictions(fit2)
#'
#' # display the average predictions
#' fit3$counterfactual.means
#'
#' @seealso [predict_counterfactuals()] for generating counterfactual
#' predictions.
#' @seealso [estimate_varcov()] for estimating the variance-covariate matrix
#' of mariginal effects
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#'
average_predictions <- function(object) {
  # assert object already stores counterfactual predictions
  if (!"counterfactual.predictions" %in% names(object)) {
    msg <- sprintf(
      "Missing counterfactual predictions. First run `%1$s <- predict_counterfactuals(%1$s)`.",
      deparse(substitute(object))
    )
    stop(msg, call. = FALSE)
  }
  trt <- attributes(object$counterfactual.predictions)$treatment.variable

  object <- .assert_sanitized(object, trt, warn = T)

  counterfactual.means <- colMeans(object$counterfactual.predictions)
  names(counterfactual.means) <- levels(.get_data(object)[[trt]])

  object$counterfactual.means <- counterfactual.means

  return(object)
}
