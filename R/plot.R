#' Plot method for beeca objects
#'
#' @description
#'
#' Creates visualizations for beeca analysis results. Currently supports
#' forest plots showing treatment effect estimates with confidence intervals.
#' Additional plot types will be added in future versions.
#'
#' @param x a beeca object (from \link[beeca]{get_marginal_effect})
#' @param type character string specifying the plot type. Currently only
#' "forest" is supported.
#' @param ... additional arguments passed to the specific plot function
#'
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   trial01$trtp <- factor(trial01$trtp)
#'   fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'     get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#'
#'   # Forest plot (default)
#'   plot(fit1)
#'
#'   # Explicit type specification
#'   plot(fit1, type = "forest")
#'
#'   # Customize
#'   plot(fit1, conf.level = 0.90, title = "My Treatment Effect")
#' }
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [plot_forest()] for forest plot details
#' @seealso [print.beeca()] for concise output
#' @seealso [summary.beeca()] for detailed summary output
#' @seealso [tidy.beeca()] for tidied parameter estimates
#' @seealso [as_gt()] for publication-ready tables
plot.beeca <- function(x, type = c("forest"), ...) {

  type <- match.arg(type)

  switch(type,
    forest = plot_forest(x, ...)
  )
}
