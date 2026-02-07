#' Forest plot for marginal treatment effects
#'
#' @importFrom rlang .data
#'
#' @description
#'
#' Creates a forest plot displaying treatment effect estimates with confidence
#' intervals. Useful for visualizing results from covariate-adjusted analyses,
#' especially with multiple treatment comparisons.
#'
#' @details
#'
#' The forest plot displays point estimates as dots with horizontal lines
#' representing confidence intervals. A vertical reference line is drawn at
#' the null effect value (0 for differences, 1 for ratios). For multiple
#' comparisons (e.g., 3-arm trials), each comparison is shown on a separate row.
#'
#' The plot can be customized using standard ggplot2 functions by adding
#' layers to the returned ggplot object.
#'
#' @param x a beeca object (from \link[beeca]{get_marginal_effect})
#' @param conf.level confidence level for confidence intervals. Default is 0.95.
#' @param title optional plot title. If NULL, a default title is generated.
#' @param xlab optional x-axis label. If NULL, a label based on the contrast
#' type is generated.
#' @param show_values logical indicating whether to display numerical values
#' on the plot. Default is TRUE.
#' @param null_line_color color for the null effect reference line.
#' Default is "darkgray".
#' @param point_size size of the point estimates. Default is 3.
#' @param ... additional arguments (not currently used)
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
#'   # Basic forest plot
#'   plot_forest(fit1)
#'
#'   # Customize with ggplot2
#'   plot_forest(fit1) +
#'     ggplot2::theme_minimal() +
#'     ggplot2::labs(title = "Treatment Effect: My Study")
#'
#'   # Risk ratio example
#'   fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'     get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
#'   plot_forest(fit_rr)
#' }
#'
#' @seealso [get_marginal_effect()] for the main analysis function
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#' @seealso [plot.beeca()] for the generic plot method
#' @seealso [tidy.beeca()] for extracting estimates in tabular form
#' @seealso [print.beeca()] for concise output
#' @seealso [summary.beeca()] for detailed summary output
#' @seealso [as_gt()] for publication-ready tables
plot_forest <- function(x,
                        conf.level = 0.95,
                        title = NULL,
                        xlab = NULL,
                        show_values = TRUE,
                        null_line_color = "darkgray",
                        point_size = 3,
                        ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for forest plots.\n",
         "Install it with: install.packages('ggplot2')",
         call. = FALSE)
  }

  # Extract tidy data
  tidy_data <- tidy(x, conf.int = TRUE, conf.level = conf.level)

  # Determine null effect value based on contrast type
  contrast_label <- attr(x$marginal_est, "contrast")
  if (is.null(contrast_label)) {
    contrast_type <- "diff"
  } else {
    contrast_type <- sub(":.*", "", contrast_label[1])
  }

  # Set null value and x-axis label based on contrast
  if (contrast_type %in% c("rr", "or")) {
    null_value <- 1
    if (is.null(xlab)) {
      xlab <- if (contrast_type == "rr") "Risk Ratio" else "Odds Ratio"
    }
  } else if (contrast_type %in% c("logrr", "logor")) {
    null_value <- 0
    if (is.null(xlab)) {
      xlab <- if (contrast_type == "logrr") "Log Risk Ratio" else "Log Odds Ratio"
    }
  } else {
    null_value <- 0
    if (is.null(xlab)) {
      xlab <- "Risk Difference"
    }
  }

  # Create default title if not provided
  if (is.null(title)) {
    title <- sprintf("Marginal Treatment Effect (%d%% CI)",
                     round(conf.level * 100))
  }

  # Reverse order for forest plot (top to bottom)
  tidy_data$term <- factor(tidy_data$term, levels = rev(tidy_data$term))

  # Create the plot
  p <- ggplot2::ggplot(tidy_data,
                       ggplot2::aes(x = .data$estimate, y = .data$term)) +
    # Null effect line
    ggplot2::geom_vline(xintercept = null_value,
                        linetype = "dashed",
                        color = null_line_color,
                        linewidth = 0.5) +
    # Confidence intervals
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$conf.low,
                                        xmax = .data$conf.high),
                           width = 0.2,
                           linewidth = 0.7,
                           orientation = "y") +
    # Point estimates
    ggplot2::geom_point(size = point_size, color = "black") +
    # Labels
    ggplot2::labs(
      title = title,
      x = xlab,
      y = "Comparison"
    ) +
    # Theme
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )

  # Add numerical values if requested
  if (show_values) {
    # Format text
    tidy_data$label_text <- sprintf("%.3f (%.3f, %.3f)",
                                     tidy_data$estimate,
                                     tidy_data$conf.low,
                                     tidy_data$conf.high)

    # Determine position (to the right of the plot)
    x_range <- range(c(tidy_data$conf.low, tidy_data$conf.high))
    x_width <- x_range[2] - x_range[1]
    text_x <- x_range[2] + x_width * 0.15

    p <- p +
      ggplot2::annotate("text",
                        x = text_x,
                        y = tidy_data$term,
                        label = tidy_data$label_text,
                        hjust = 0,
                        size = 3.5)

    # Extend x-axis to make room for text
    p <- p +
      ggplot2::coord_cartesian(xlim = c(x_range[1],
                                        x_range[2] + x_width * 0.4),
                               clip = "off")
  }

  p
}
