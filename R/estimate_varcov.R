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
#' For more details, see Magirr et al. (2025) \doi{10.1002/pst.70021}.
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
#' @seealso [predict_counterfactuals()] for generating counterfactual predictions
#' @seealso [average_predictions()] for averaging counterfactual predictions
#' @seealso [apply_contrast()] for computing a summary measure
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#' @seealso [beeca_fit()] for streamlined analysis pipeline
#'
#' @references Ye, T., Shao, J., Yi, Y., and Zhao, Q. (2023). "Robust Variance
#' Estimation for Covariate-Adjusted Unconditional Treatment Effect in Randomized
#' Clinical Trials with Binary Outcomes." *Statistical Theory and Related Fields*
#' 7(2): 159-163. <https://doi.org/10.1080/24754269.2023.2205802>
#'
#' @references Ge, M., Durham, L. K., Meyer, R. D., Xie, W., and Flournoy, N. (2011).
#' "Covariate-Adjusted Difference in Proportions from Clinical Trials Using Logistic
#' Regression and Weighted Risk Differences." *Drug Information Journal* 45(4): 481-493.
#' <https://doi.org/10.1177/009286151104500409>
#'
#' @references Bannick, M. S., Jun, Y., Josey, K., Weinberg, C. R., and Karlson, E. W.
#' (2023). "A General Form of Covariate Adjustment in Randomized Clinical Trials."
#' arXiv preprint arXiv:2306.10213. <https://arxiv.org/abs/2306.10213>
#'
estimate_varcov <- function(object, strata = NULL, method = c("Ge", "Ye"),
                            type = c("HC0", "model-based", "HC3", "HC", "HC1", "HC2", "HC4", "HC4m", "HC5"),
                            mod = FALSE) {
  # Assert means are available in object
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
    # Throw warning is strata supplied but method == "Ge" is specified
    if (!is.null(strata)) {
      warning("Strata is supplied but will be ignored when using method='Ge'", call. = FALSE)
    }
    varcov <- varcov_ge(object, trt, type)
  }

  # Assert varcov is symmetric
  if (isFALSE(all.equal(varcov[lower.tri(varcov)], varcov[upper.tri(varcov)]))) {
    stop("Error in calculation of variance covariance matrix: `varcov` is not symmetric.",
         call. = FALSE)
  }

  object$robust_varcov <- varcov
  attr(object$robust_varcov, "type") <- ifelse(method == "Ge", paste0(method, " - ", type), method)

  return(object)
}


## Re-implementation based on Ye et al. (https://doi.org/10.1080/24754269.2023.2205802)
varcov_ye <- function(object, trt, strata, mod) {

  data <- .get_data(object)

  # Extract allocation information
  pi_all <- prop.table(table(data[, trt]))
  n <- length(data[, trt])

  # Extract observed responses
  y_observed <- lapply(levels(data[[trt]]), \(x) object$y[data[, trt] == x])
  names(y_observed) <- levels(data[[trt]])

  # Extract counterfactual predictions
  cf_pred <- object$counterfactual.predictions
  # Append actual treatment assignment to cf predictions
  cf_pred[[trt]] <- data[,trt]

  # Estimate varcov using Ye et al (2022)

  # diag_vars stores variances per treatment arm
  diag_vars <- list()
  if (mod) {

    for (trtlvl in levels(data[[trt]])) {
      var_resid <- var(y_observed[[trtlvl]] - cf_pred[cf_pred[[trt]] == trtlvl, trtlvl])
      cov_obs_pred <- cov(y_observed[[trtlvl]], cf_pred[cf_pred[[trt]] == trtlvl, trtlvl])
      var_pred <- var(cf_pred[, trtlvl])

      diag_vars[[trtlvl]] <- var_resid / pi_all[[trtlvl]] + 2 * cov_obs_pred - var_pred

    }

  } else {
    # Use variance decomposition to match RobinCar implementation

    for (trtlvl in levels(data[[trt]])) {
      var_obs <- var(y_observed[[trtlvl]])
      cov_obs_pred <- cov(y_observed[[trtlvl]], cf_pred[cf_pred[[trt]] == trtlvl, trtlvl])
      var_pred <- var(cf_pred[, trtlvl])

      var_obs_pred <- var_obs - 2 * cov_obs_pred + var_pred

      diag_vars[[trtlvl]] <- var_obs_pred / pi_all[[trtlvl]] + 2 * cov_obs_pred - var_pred

    }

  }

  # tri_vars stores trt arm covariances
  tri_vars <- list()
  for (trtlvl in levels(data[[trt]])) {

    # Covar of observed arm i and corresponding counterfactual predictions for those on arm j
    cov_obsi_predj <- sapply(levels(data[[trt]])[levels(data[[trt]]) != trtlvl],
                             \(x) cov(y_observed[[trtlvl]], cf_pred[cf_pred[[trt]] == trtlvl, x]))

    # Covar of observed arm j and corresponding counterfactual predictions for those on arm i
    cov_obsj_predi <- sapply(levels(data[[trt]])[levels(data[[trt]]) != trtlvl],
                             \(x) cov(y_observed[[x]], cf_pred[cf_pred[[trt]] == x, trtlvl]))

    # Covar of counterfactual predictions for arm i and arm j for all subjects
    cov_predi_predj <- sapply(levels(data[[trt]])[levels(data[[trt]]) != trtlvl],
                              \(x) cov(cf_pred[, trtlvl], cf_pred[, x]))

    tri_vars[[trtlvl]] <- cov_obsi_predj + cov_obsj_predi - cov_predi_predj

  }

  # Construct variance covariance matrix
  varcov <- matrix(nrow = nlevels(data[[trt]]),
                   ncol = nlevels(data[[trt]]))
  rownames(varcov) <- colnames(varcov) <- levels(data[[trt]])

  for (trtlvl in rownames(varcov)) {
    varcov[trtlvl, c(trtlvl, names(tri_vars[[trtlvl]]))] <- c(diag_vars[[trtlvl]], tri_vars[[trtlvl]])
  }

  # Variance adjustment based on stratified CAR scheme, formula 18: https://doi.org/10.1080/01621459.2022.2049278
  if (!is.null(strata)) {

    # Warn if any single stratification variable has suspiciously many unique values
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

  data <- .get_data(object)

  # Assert strata is a column name found in model data
  if (!all(strata %in% colnames(data))) {
    msg <- sprintf('Stratification variable(s) "%s" must be found in the model', paste(strata[!strata %in% colnames(data)], collapse = ", "))
    stop(msg, call. = FALSE)
  }

  strata_all <- Reduce(function(...) paste(..., sep = "-"), data[strata])
  data <- cbind(data, strata_all)

  # Extract allocation information
  pi_all <- prop.table(table(data[,trt]))
  # Get predictions for observed data and calculate residuals
  data$preds <- predict(object, type = "response")
  data$resid <- object$y - data$preds

  get_rb <- function(z) {
    rb_list <- list()
    for (trtlvl in levels(data[[trt]])) {
      resid_i_z <- data$resid[(data[[trt]] == trtlvl) & (data$strata_all == z)]
      rb_list[[trtlvl]] <- (1 / pi_all[[trtlvl]]) * mean(resid_i_z)
    }
    return(diag(rb_list))
  }

  omega_sr <- diag(pi_all) - pi_all %*% t(pi_all)

  # Always assume permuted block stratified randomization
  omega_z <- matrix(0, nrow = nlevels(data[[trt]]), ncol = nlevels(data[[trt]]))

  # Compute outer expectation
  strata_levels <- levels(as.factor(data$strata_all))
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

  # Counterfactual predictions
  cf_pred <- object$counterfactual.predictions

  # Get derivatives per treatment level counterfactual outcomes
  d_list <- list()
  for (trtlvl in levels(data[[trt]])) {
    X_i <- data
    X_i[, trt] <- factor(trtlvl, levels = levels(data[[trt]]))
    X_i <- model.matrix(object$formula, X_i)

    # Compute derivative
    pderiv_i <- cf_pred[, trtlvl] * (1 - cf_pred[, trtlvl])
    d_i <- (t(pderiv_i) %*% as.matrix(X_i)) / nrow(X_i)
    d_list[[trtlvl]] <- d_i
  }

  all_d <- do.call(rbind, d_list)
  varcov_ge <- all_d %*% V %*% t(all_d)
  rownames(varcov_ge) <- colnames(varcov_ge) <- levels(data[[trt]])

  return(varcov_ge)
}
