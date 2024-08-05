#' Augment original data with counter factual predictions
#'
#'
#' @examples
#' ## Set treatment to a factor
#' trial01$trtp <- factor(trial01$trtp)
#'
#' ## Fit a logistic regression working model and pass it to beeca
#' fit1 <- glm(aval ~ trtp + bl_cov, family="binomial", data=trial01) |>
#'   get_marginal_effect(trt="trtp", method="Ye", contrast="diff")
#'
#'   augment_beeca(fit1)
#'
#' @export
augment_beeca <- function(x) {

  # Get the data and predictions and combine
  model_data <- .get_data(x)
  predictions <- x$counterfactual.predictions
  data <- dplyr::bind_cols(model_data, predictions)

  # merge with original data in case patients dropped due to missinginess
  original <- x$data
  data <- original |>
    dplyr::left_join(data) |>
    dplyr::as_tibble()

  # todo: mutate additional statistics for CI table - check text book

  return(data)
}

