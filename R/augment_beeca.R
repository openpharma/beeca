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

