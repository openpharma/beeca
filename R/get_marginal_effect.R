#' Estimate marginal treatment effects using a GLM working model
#'
#' @description
#'
#' Estimates the marginal treatment effect from a logistic regression working model
#' using a specified choice of variance estimator and contrast.
#'
#' @details
#' The `get_marginal_effect` function is a wrapper that facilitates
#' advanced variance estimation techniques for GLM models with covariate adjustment
#' targeting a population average treatment effect. It is particularly useful in clinical trial analysis
#' and other fields requiring robust statistical inference.
#' It allows researchers to account for complex study designs,
#' including stratification and treatment contrasts, by providing a flexible
#' interface for variance-covariance estimation.
#'
#' @param object a fitted \link[stats]{glm} object.
#' @param trt a string specifying the name of the treatment variable
#' in the model formula. It must be one of the linear predictor variables used
#' in fitting the `object`.
#' @param strata an optional string or vector of strings specifying the names
#' of stratification variables. Relevant only for Ye's method and used to
#' adjust the variance-covariance estimation for stratification.
#' If provided, each specified variable must be present in the model.
#' @param method a string indicating the chosen method for variance estimation.
#' Supported methods are `Ge` and `Ye`. The default method is `Ge` based on Ge
#' et al (2011) which is suitable for the variance estimation of conditional
#' average treatment effect. The method `Ye` is based on Ye et al (2023) and is
#' suitable for the variance estimation of population average treatment effect.
#' For more details, see Magirr et al. (2025) \doi{10.1002/pst.70021}.
#'
#' @param type a string indicating the type of
#' variance estimator to use (only applicable for Ge's method). Supported types include HC0 (default),
#' model-based, HC3, HC, HC1, HC2, HC4, HC4m, and HC5. See \link[sandwich]{vcovHC} for heteroscedasticity-consistent estimators.
#'
#' @param contrast a string indicating choice of contrast. Defaults to 'diff' for a risk difference. See \link[beeca]{apply_contrast}.
#'
#' @param reference a string or list of strings indicating which treatment
#' group(s) to use as reference level for pairwise comparisons. Accepted values
#' must be a subset of the levels in the treatment variable. Default to the
#' first n-1 treatment levels used in the `glm` object. This parameter influences
#' the calculation of treatment effects relative to the chosen reference group.
#'
#' @param mod for Ye's method, the implementation of open-source RobinCar package
#' has an additional variance decomposition step when estimating the robust variance,
#' which then utilizes different counterfactual outcomes than the original reference.
#' Set `mod = TRUE` to use exactly the implementation method described in
#' Ye et al (2022), default to `FALSE` to use the modified implementation in
#' RobinCar and Bannick et al (2023) which improves stability.
#'
#' @return an updated `glm` object appended with marginal estimate components:
#' counterfactual.predictions (see \link[beeca]{predict_counterfactuals}),
#' counterfactual.means (see \link[beeca]{average_predictions}),
#' robust_varcov (see \link[beeca]{estimate_varcov}),
#' marginal_est, marginal_se (see \link[beeca]{apply_contrast}) and marginal_results. A summary is shown below
#' \tabular{ll}{
#'  counterfactual.predictions \tab Counterfactual predictions based on the working model. For each subject in the input glm data, the potential outcomes are obtained by assigning subjects to each of the possible treatment variable levels. Each prediction is associated with a descriptive label explaining the counterfactual scenario. \cr
#'  counterfactual.means       \tab Average of the counterfactual predictions for each level of the treatment variable. \cr
#'  robust_varcov              \tab Variance-covariance matrix of the marginal effect estimate for each level of treatment variable, with estimation method indicated in the attributes.\cr
#'  marginal_est               \tab Marginal treatment effect estimate for a given contrast.\cr
#'  marginal_se                \tab Standard error estimate of the marginal treatment effect estimate. \cr
#'  marginal_results           \tab Analysis results data (ARD) containing a summary of the analysis for subsequent reporting. \cr
#' }
#' @seealso [predict_counterfactuals()] for generating counterfactual predictions
#' @seealso [average_predictions()] for averaging counterfactual predictions
#' @seealso [estimate_varcov()] for robust variance estimation
#' @seealso [apply_contrast()] for computing treatment contrasts
#' @seealso [beeca_fit()] for streamlined convenience wrapper
#' @seealso [tidy.beeca()] for tidied parameter estimates
#' @seealso [summary.beeca()] for detailed summary output
#' @seealso [print.beeca()] for concise output
#' @seealso [plot.beeca()] and [plot_forest()] for visualizations
#' @seealso [augment.beeca()] for augmented data with predictions
#' @seealso [as_gt()] for publication-ready tables
#' @seealso [beeca_to_cards_ard()] for cards ARD integration
#'
#' @importFrom utils packageVersion
#' @export
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
#' fit1$marginal_results
get_marginal_effect <- function(object, trt, strata = NULL,
                                method = "Ge",
                                type = "HC0",
                                contrast = "diff",
                                reference, mod = FALSE) {
  object <- .assert_sanitized(object, trt)

  object <- object |>
    predict_counterfactuals(trt) |>
    average_predictions() |>
    estimate_varcov(strata, method, type, mod) |>
    apply_contrast(contrast, reference)

  data <- .get_data(object)
  n_trt_levels <- nlevels(data[[trt]])
  outcome <- paste(object$formula[[2]])

  # Marginal results for each trt level
  marginal_responses <- data.frame(
    TRTVAR = rep(trt, 5*n_trt_levels),
    TRTVAL = rep(levels(data[[trt]]), each=5),
    PARAM = rep(outcome, 5),
    ANALTYP1 = rep(c(rep("DESCRIPTIVE", 3), rep("INFERENTIAL", 2)), n_trt_levels),
    STAT = rep(c("N", "n", "%", "risk", "risk_se"), n_trt_levels),
    STATVAL = c(sapply(levels(data[[trt]]), \(x) {
      c(
        sum(data[trt] == x),
        sum(data[outcome][data[trt] == x] == "1"),
        sum(data[outcome][data[trt] == x] == "1") / sum(data[trt] == x) * 100,
        object$counterfactual.means[x],
        sqrt(object$robust_varcov[x, x])
      )
    })
    ),
    ANALMETH = rep(c("count", "count", "percentage", "g-computation",
                     attributes(object$marginal_se)[["type"]]), n_trt_levels)
  )

  # Contrast results
  n_contrasts <- length(object$marginal_est)

  marginal_contrasts <- data.frame(TRTVAR = rep(trt, 2*n_contrasts),
                                   TRTVAL = c(rbind(attributes(object$marginal_est)[["contrast"]],
                                                    attributes(object$marginal_se)[["contrast"]])),
                                   PARAM = rep(outcome, 2*n_contrasts),
                                   ANALTYP1 = c(rep("INFERENTIAL", 2*n_contrasts)),
                                   STAT = rep(c(contrast, paste0(contrast, "_se")), n_contrasts),
                                   STATVAL = c(rbind(object$marginal_est, object$marginal_se)),
                                   ANALMETH = rep(c("g-computation", attributes(object$marginal_se)[["type"]]), n_contrasts)
  )

  object$marginal_results <- rbind(marginal_responses,
                                   marginal_contrasts) |> dplyr::as_tibble()
  object$marginal_results$ANALDESC <- paste0("Computed using beeca@", packageVersion("beeca"))

  # Add beeca class to enable S3 methods (tidy, augment, etc.)
  class(object) <- c("beeca", class(object))

  return(object)
}
