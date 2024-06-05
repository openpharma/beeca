#' Estimate variance-covariance matrix for marginal estimand based on GLM model
#'
#' @description
#'
#' Main variance estimation function. Estimates the variance-covariance
#' matrix of a marginal estimand for a generalized linear model (GLM) object
#' using specified methods. This function supports both Ge's and Ye's methods
#' for variance estimation, accommodating different estimand specifications.
#' @importFrom stats cov var
#' @param object a fitted \code{\link[stats]{glm}} object augmented with
#' `counterfactual.predictions`, `counterfactual.predictions` and `counterfactual.means`
#'
#' @param strata an optional string or vector of strings specifying the names
#' of stratification variables. Relevant only for Ye's method and used to
#' adjust the variance-covariance estimation for stratification.
#' If provided, each specified variable must be present in the model.
#'
#' @param method a string indicating the chosen method for variance estimation.
#' Supported methods are `Ge` and `Ye`. The default method is `Ge` based on Ge
#' et al (2011) which is suitable for the variance estimation of conditional
#' average treatment effect. The method `Ye` is based on Ye et al (2023) and is
#' suitable for the variance estimation of population average treatment effect.
#' For more details, see [Magirr et al. (2024)](https://osf.io/9mp58/).
#'
#'
#' @param type a string indicating the type of
#' variance estimator to use (only applicable for Ge's method). Supported types include HC0 (default),
#' model-based, HC3, HC, HC1, HC2, HC4, HC4m, and HC5. See \link[sandwich]{vcovHC} for heteroscedasticity-consistent estimators.
#' This parameter allows for flexibility in handling heteroscedasticity
#' and model specification errors.
#'
#' @param mod For Ye's method, the implementation of open-source RobinCar package
#' has an additional variance decomposition step when estimating the robust variance,
#' which then utilizes different counterfactual outcomes than the original reference.
#' Set `mod = TRUE` to use exactly the implementation method described in
#' Ye et al (2022), default to `FALSE` to use the modified implementation in
#' RobinCar and Bannick et al (2023) which improves stability.
#'
#' @return an updated `glm` object appended with an
#' additional component `robust_varcov`, which is the estimated variance-covariance matrix
#'  of the marginal effect. The matrix format and estimation method are
#'  indicated in the matrix attributes.
#'
#' @details
#' The `estimate_varcov` function facilitates robust variance estimation
#' techniques for GLM models, particularly useful in clinical trial analysis
#' and other fields requiring robust statistical inference. It allows
#' researchers to account for complex study designs,
#' including stratification and different treatment contrasts,
#' by providing a flexible interface for variance-covariance estimation.
#'
#' Note: Ensure that the `glm` object has been adequately prepared with
#' \code{\link{predict_counterfactuals}} and \code{\link{average_predictions}}
#' before applying `estimate_varcov()`. Failure to do so may result in
#' errors indicating missing components.
#'
#' @examples
#' # Example usage with a binary outcome GLM model
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
#'
#' #' # Preprocess fit1 as required by estimate_varcov
#' fit2 <- fit1 |>
#'   predict_counterfactuals(trt = "trtp") |>
#'   average_predictions()
#'
#' # Estimate variance-covariance using Ge's method
#' fit3_ge <- estimate_varcov(fit2, method = "Ge")
#' print(fit3_ge$robust_varcov)
#'
#'
#' # Estimate variance-covariance using Ye's method with stratification
#' fit4 <- glm(aval ~ trtp + bl_cov_c, family = "binomial", data = trial01) |>
#'   predict_counterfactuals(trt = "trtp") |>
#'   average_predictions()
#' fit4_ye <- estimate_varcov(fit4, method = "Ye", strata = "bl_cov_c")
#' print(fit4_ye$robust_varcov)
#'
#' @export
#'
#' @seealso [average_predictions()] for averaging counterfactual
#' predictions.
#' @seealso [apply_contrast()] for computing a summary measure.
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#'
#' @references Ye T. et al. (2023) Robust variance
#' estimation for covariate-adjusted unconditional treatment effect in randomized
#' clinical trials with binary outcomes. Statistical Theory and Related Fields
#'
#' @references Ge M. et al. (2011) Covariate-Adjusted
#' Difference in Proportions from Clinical Trials Using Logistic Regression
#' and Weighted Risk Differences. Drug Information Journal.
#'
#' @references Bannick, M. S., et al. A General Form of Covariate Adjustment in
#' Randomized Clinical Trials. arXiv preprint arXiv:2306.10213 (2023).
#'
estimate_varcov <- function(object, strata = NULL, method = c("Ge", "Ye"),
                            type = c("HC0", "model-based", "HC3", "HC", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                            mod = FALSE) {
  # assert means are available in object
  if (!"counterfactual.means" %in% names(object)) {
    msg <- sprintf(
      "Missing counterfactual means. First run `%1$s <- average_predictions(%1$s)`.",
      deparse(substitute(object))
    )
    stop(msg, call. = FALSE)
  }

  trt <- attributes(object$counterfactual.predictions)$treatment.variable
  object <- .assert_sanitized(object, trt, warn = T)

  method <- match.arg(method)
  if (method == "Ye") {
    varcov <- varcov_ye(object, trt, strata, mod)
  } else if (method == "Ge") {
    type <- match.arg(type)
    # throw warning is strata supplied but method == "Ge" is specified
    if (!is.null(strata)) {
      warning("Strata is supplied but will be ignored when using method='Ge'", call. = FALSE)
    }
    varcov <- varcov_ge(object, trt, type)
  }

  rownames(varcov) <- colnames(varcov) <- levels(.get_data(object)[[trt]])

  object$robust_varcov <- varcov
  attr(object$robust_varcov, "type") <- ifelse(method == "Ge", paste0(method, " - ", type), method)

  return(object)
}


## Re-implementation based on Ye et al. (https://doi.org/10.1080/24754269.2023.2205802)
varcov_ye <- function(object, trt, strata, mod) {
  ## calculate quantities
  # get data
  data <- .get_data(object)

  # counterfactual predictions
  cf_pred <- object$counterfactual.predictions

  pred_0 <- cf_pred$cf_pred_0
  pred_1 <- cf_pred$cf_pred_1

  ## Extract allocation information
  pi_0 <- mean(data[, trt] == levels(data[[trt]])[1])
  pi_1 <- 1 - pi_0
  n <- length(data[, trt])

  ## Extract responses on observed arm i and arm j
  y_0 <- object$y[data[, trt] == levels(data[[trt]])[1]]
  y_1 <- object$y[data[, trt] == levels(data[[trt]])[2]]

  ## Predictions on arm i from patients in arm i only
  pred_0_0 <- pred_0[data[, trt] == levels(data[[trt]])[1]]
  pred_1_1 <- pred_1[data[, trt] == levels(data[[trt]])[2]]

  ## Predictions on arm i from patients in arm j only
  pred_0_1 <- pred_0[data[, trt] == levels(data[[trt]])[2]]
  pred_1_0 <- pred_1[data[, trt] == levels(data[[trt]])[1]]

  ## estimate varcov using Ye et al (2022)
  if (mod) {
    v_0_0 <- var(y_0 - pred_0_0) / pi_0 + 2 * cov(y_0, pred_0_0) - var(pred_0)
    v_1_1 <- var(y_1 - pred_1_1) / pi_1 + 2 * cov(y_1, pred_1_1) - var(pred_1)
  } else {
    # use variance decomposition to match RobinCar implementation
    var_y0_pred00 <- var(y_0) - 2 * cov(y_0, pred_0_0) + var(pred_0)
    var_y1_pred11 <- var(y_1) - 2 * cov(y_1, pred_1_1) + var(pred_1)

    v_0_0 <- var_y0_pred00 / pi_0 + 2 * cov(y_0, pred_0_0) - var(pred_0)
    v_1_1 <- var_y1_pred11 / pi_1 + 2 * cov(y_1, pred_1_1) - var(pred_1)
  }

  v_0_1 <- v_1_0 <- cov(y_0, pred_1_0) + cov(y_1, pred_0_1) - cov(pred_0, pred_1)

  ## Variance of theta_i = p_i(Y=1) on arm i = 0,1
  varcov <- matrix(c(
    v_0_0, v_0_1,
    v_1_0, v_1_1
  ), nrow = 2, ncol = 2)

  ## Variance adjustment based on stratified CAR scheme, formula 18: https://doi.org/10.1080/01621459.2022.2049278
  if (!is.null(strata)) {

    # warn if any single stratification variable has suspiciously many unique values
    n.strata <- sapply(strata, function(x) length(unique(.get_data(object)[[x]])))
    if (any(n.strata > 3)) {
      msg <- sprintf(
        "More than three unique values are found in stratification variable %s. Please double check if the correct stratification variables are supplied via `strata` argument.",
        paste(strata[which(n.strata > 3)], collapse = ", ")
      )
      warning(msg, call. = FALSE)
    }

    varcov <- varcov - ye_strata_adjustment(object, trt, strata)
  }
  varcov <- varcov / n
  return(varcov)
}

ye_strata_adjustment <- function(object, trt, strata) {
  # get data
  data <- .get_data(object)

  # assert strata is a column name found in model data
  if (!all(strata %in% colnames(data))) {
    msg <- sprintf('Stratification variable(s) "%s" must be found in the model', paste(strata[!strata %in% colnames(data)], collapse = ", "))
    stop(msg, call. = FALSE)
  }

  strata_all <- Reduce(function(...) paste(..., sep = "-"), data[strata])
  data <- cbind(data, strata_all)

  pi_0 <- mean(data[, trt] == levels(data[[trt]])[1])
  pi_1 <- 1 - pi_0
  data$preds <- predict(object, type = "response")
  data$resid <- object$y - data$preds

  get_rb <- function(z) {
    # assume always two treatment levels
    resid_0_z <- data$resid[(data[, trt] == levels(data[[trt]])[1]) & (data$strata_all == z)]
    resid_1_z <- data$resid[(data[, trt] == levels(data[[trt]])[2]) & (data$strata_all == z)]
    rb_0 <- (1 / pi_0) * mean(resid_0_z)
    rb_1 <- (1 / pi_1) * mean(resid_1_z)
    return(diag(c(rb_0, rb_1)))
  }

  strata_levels <- levels(as.factor(data$strata_all))

  omega_sr <- diag(c(pi_0, pi_1)) - c(pi_0, pi_1) %*% t(c(pi_0, pi_1))
  # always assume permuted block stratified randomization
  omega_z <- matrix(0, nrow = length(unique((data[[trt]]))), ncol = length(unique((data[[trt]]))))

  # compute outer expectation
  erb_s <- lapply(strata_levels, function(x) {
    erb_i <- get_rb(x) %*% (omega_sr - omega_z) %*% get_rb(x)
    erb_i * (sum(data$strata_all == x) / nrow(data))
  })
  erb <- Reduce("+", erb_s)

  return(erb)
}




## Re-implementation of Ge et al. (https://doi.org/10.1177/009286151104500409)
#' @importFrom stats vcov
varcov_ge <- function(object, trt, type) {
  data <- .get_data(object)

  # Sandwich estimator of variance-covariance matrix
  if (type == "model-based") {
    V <- stats::vcov(object)
  } else {
    V <- sandwich::vcovHC(object, type = type)
  }

  # counterfactual predictions
  cf_pred <- object$counterfactual.predictions

  pred_0 <- cf_pred$cf_pred_0
  pred_1 <- cf_pred$cf_pred_1

  # get order of treatment column and assign counterfactual treatments
  X0 <- X1 <- data
  X0[, trt] <- factor(levels(data[[trt]])[1], levels = levels(data[[trt]]))
  X1[, trt] <- factor(levels(data[[trt]])[2], levels = levels(data[[trt]]))
  X0 <- model.matrix(object$formula, X0)
  X1 <- model.matrix(object$formula, X1)
  # compute derivatives
  pderiv0 <- pred_0 * (1 - pred_0)
  pderiv1 <- pred_1 * (1 - pred_1)
  d0 <- (t(pderiv0) %*% as.matrix(X0)) / nrow(X0)
  d1 <- (t(pderiv1) %*% as.matrix(X1)) / nrow(X1)

  # variance components aligned with varcov_ye implementation
  v_0_0 <- d0 %*% V %*% t(d0)
  v_1_1 <- d1 %*% V %*% t(d1)
  v_0_1 <- d0 %*% V %*% t(d1)
  v_1_0 <- d1 %*% V %*% t(d0)

  varcov_ge <- matrix(
    c(
      v_0_0, v_0_1,
      v_1_0, v_1_1
    ),
    nrow = 2, ncol = 2
  )
  return(varcov_ge)
}
