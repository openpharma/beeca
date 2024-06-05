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
#' For more details, see [Magirr et al. (2024)](https://osf.io/9mp58/).
#'
#' @param type a string indicating the type of
#' variance estimator to use (only applicable for Ge's method). Supported types include HC0 (default),
#' model-based, HC3, HC, HC1, HC2, HC4, HC4m, and HC5. See \link[sandwich]{vcovHC} for heteroscedasticity-consistent estimators.
#'
#' @param contrast a string indicating choice of contrast. Defaults to 'diff' for a risk difference. See \link[beeca]{apply_contrast}.
#'
#' @param reference a string indicating which treatment group should be considered as
#' the reference level. Accepted values are one of the levels in the treatment
#' variable. Default to the first level used in the `glm` object. This parameter influences the calculation of treatment effects
#' relative to the chosen reference group.
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

  trt_ref <- attr(object$marginal_est, "reference")
  trt_inv <- levels(data[[trt]])[-which(levels(data[[trt]]) == trt_ref)]
  outcome <- paste(object$formula[[2]])

  marginal_results <- data.frame(
    TRTVAR = c(rep(trt, 12)),
    TRTVAL = c(
      rep(trt_inv, 5), rep(trt_ref, 5),
      rep(attributes(object$marginal_est)[["contrast"]], 2)
    ),
    PARAM = rep(outcome, 12),
    ANALTYP1 = c(rep(c(rep("DESCRIPTIVE", 3), rep("INFERENTIAL", 2)), 2), rep("INFERENTIAL", 2)),
    STAT = c("N", "n", "%", "risk", "risk_se", "N", "n", "%", "risk", "risk_se", contrast, paste0(contrast, "_se")),
    STATVAL = c(
      sum(data[trt] == trt_inv), sum(data[outcome][data[trt] == trt_inv] == "1"),
      sum(data[outcome][data[trt] == trt_inv] == "1") / sum(data[trt] == trt_inv) * 100,
      object$counterfactual.means[trt_inv], sqrt(object$robust_varcov[trt_inv, trt_inv]),
      sum(data[trt] == trt_ref), sum(data[outcome][data[trt] == trt_ref] == "1"),
      sum(data[outcome][data[trt] == trt_ref] == "1") / sum(data[trt] == trt_ref) * 100,
      object$counterfactual.means[trt_ref], sqrt(object$robust_varcov[trt_ref, trt_ref]),
      object$marginal_est, object$marginal_se
    ),
    ANALMETH = c(
      rep(c("count", "count", "percentage", "g-computation", attributes(object$marginal_se)[["type"]]), 2),
      "g-computation", attributes(object$marginal_se)[["type"]]
    ),
    ANALDESC = paste0("Computed using beeca@", packageVersion("beeca"))
  ) |>
    dplyr::as_tibble()

  object$marginal_results <- marginal_results

  return(object)
}
