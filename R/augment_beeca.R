#' @export
augment_beeca <- function(x) {
  data <- .get_data(x)

  data <- dplyr::as_tibble(data)

  return(data)
}

