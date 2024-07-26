#' Tidy a beeca object
#'
#' The broom package exports a tidier for `"glm"` objects.
#' This function adds on top of that and returns more information
#' that is stored from gcomputation.
#'
#' The function also utilizes additional information stored when the
#' glm object is created with `get_marginal_effect()`.
#' It's recommended to always use this function with `get_marginal_effect()`.
#'
#' @param x a 'glm' object created with `get_marginal_effect()`
#' @param conf.int a logical indicating whether to include confidence intervals
#' @param conf.level the confidence level to use for the confidence interval
#'
#' @return a tibble
#' @export
#'
#' @examples
#' # example use of tidy_beeca
#' # Set treatment to a factor
#' trial01$trtp <- factor(trial01$trtp)
#'
#' ## Fit a logistic regression working model and pass it to beeca
#' fit1 <- glm(aval ~ trtp + bl_cov, family="binomial", data=trial01) |>
#'   get_marginal_effect(trt="trtp", method="Ye", contrast="diff")
#'
#' # compare with broom tidy. one works on logistic model(working model)
#' # tidy_beeca tidies the analyses from g-computation
#' broom::tidy(fit1)
#' tidy_beeca(fit1)
#'
#' # with confidence intervals
#' broom::tidy(fit1, conf.int = TRUE)
#' tidy_beeca(fit1, conf.int = TRUE)
#'
tidy_beeca <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {

  results <- NULL

  # check inputs ---------------------------------------------------------------

  # check if x is a beeca object
  if (!inherits(x, "beeca")) {
    stop("x must be a 'beeca' object")
  }

  results <- broom::tidy(x)

  ## Extract results
  marginal_results <- x$marginal_results

  # define a tidy tibble
  result <- tibble::tibble(
    term = marginal_results$TRTVAR[1],
    contrast = marginal_results[marginal_results$STAT == "diff", "TRTVAL"][[1]],
    estimate = marginal_results[marginal_results$STAT == "diff", "STATVAL"][[1]],
    std.error = marginal_results[marginal_results$STAT == "diff_se", "STATVAL"][[1]],
    statistic = NA,
    p.value = NA
  )

  ## add confidence interval if specified
  if (conf.int) {

    result <- result |>
      dplyr::mutate(
        conf.low = estimate - qnorm(1 - conf.level / 2) * std.error,
        conf.high = estimate + qnorm(1 - conf.level / 2) * std.error
      )
  }

  # return tidy tibble
  return(result)
}
