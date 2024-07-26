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
  data <- .get_data(x)

  data <- dplyr::as_tibble(data)
  predictions <- dplyr::as_tibble(x$counterfactual.predictions)
  data <- dplyr::bind_cols(data, predictions)

  original <- x$data

  data <- original |> dplyr::left_join(data)


  return(data)
}

