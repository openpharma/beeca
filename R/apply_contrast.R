#' Apply contrast to calculate marginal estimate of treatment effect and corresponding standard error
#'
#' @description
#'
#' Calculates the marginal estimate of treatment effect and its corresponding
#' standard error based on a fitted GLM object using specified contrast (summary
#' measure) methods
#'
#' @details
#' The `apply_constrast()` functions computes the summary measure between two arms
#' based on the estimated marginal effect and its variance-covariance matrix using
#' the Delta method.
#'
#' Note: Ensure that the `glm` object has been adequately prepared with
#' `average_predictions()` and `estimate_varcov()`
#' before applying `apply_contrast()`. Failure to do so may result in
#' errors indicating missing components.
#'
#' @param object a fitted \code{\link[stats]{glm}} object augmented with
#' `counterfactual.predictions`, `counterfactual.means` and `robust_varcov`.
#'
#' @param contrast a string specifying the type of contrast to apply.
#' Accepted values are "diff" (risk difference), "rr" (risk ratio),
#' "or" (odds ratio), "logrr" (log risk ratio), "logor" (log odds ratio).
#' Note: log-transformed ratios (logrr and logor) work better compared to rr
#' and or when computing confidence intervals using normal approximation.
#' The choice of contrast affects how treatment effects are calculated and
#' interpreted. Default is `diff`.
#'
#' @param reference a string indicating which treatment group should be considered as
#' the reference level. Accepted values are one of the levels in the treatment
#' variable. Default to the first level used in the `glm` object.
#'
#' This parameter influences the calculation of treatment effects
#' relative to the chosen reference group.
#'
#' @return An updated `glm` object with two additional components
#' appended: `marginal_est` (marginal estimate of the treatment effect)
#' and `marginal_se` (standard error of the marginal estimate).
#' These appended component provide crucial information for interpreting
#' the treatment effect using the specified contrast method.
#'
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#'
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
#'   predict_counterfactuals(trt = "trtp") |>
#'   average_predictions() |>
#'   estimate_varcov(method = "Ye") |>
#'   apply_contrast("diff", reference = "0")
#'
#' # Assuming `trial01` is a dataset with treatment (`trtp`)
#' # and baseline covariate (`bl_cov`)
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
#'
#' # Preprocess fit1 as required by apply_contrast
#' fit2 <- fit1 |>
#'   predict_counterfactuals(trt = "trtp") |>
#'   average_predictions() |>
#'   estimate_varcov(method = "Ye")
#'
#' # Apply contrast to calculate marginal estimates
#' fit3 <- apply_contrast(fit2, contrast = "diff", reference = "0")
#'
#' fit3$marginal_est
#' fit3$marginal_se
#'
#' @export
#'
apply_contrast <- function(object, contrast = c("diff", "rr", "or", "logrr", "logor"), reference) {
  # assert means are available in object
  if (!"counterfactual.means" %in% names(object)) {
    msg <- sprintf(
      "Missing counterfactual means. First run `%1$s <- average_predictions(%1$s)`.",
      deparse(quote(object))
    )
    stop(msg, call. = FALSE)
  }
  # assert varcov is available in object
  if (!"robust_varcov" %in% names(object)) {
    msg <- sprintf(
      "Missing robust varcov. First run `%1$s <- get_varcov(%1$s, ...)`.",
      deparse(quote(object))
    )
    stop(msg, call. = FALSE)
  }

  trt <- attributes(object$counterfactual.predictions)$treatment.variable
  object <- .assert_sanitized(object, trt, warn = T)

  data <- .get_data(object)
  # assert non-missing reference
  if (missing(reference)) {
    reference <- object$xlevels[[trt]][1]
    warning(sprintf("No reference argument was provided, using %s as the reference level", reference), call. = FALSE)
  }

  # assert reference to be one of the treatment levels
  if (!reference %in% levels(data[[trt]])) {
    stop("Reference must be one of : ", paste(levels(data[[trt]]), collapse = ", "), ".")
  }

  # match contrast argument and retrieve relevant functions
  contrast <- match.arg(contrast)
  c_est <- get(contrast)
  c_str <- get(paste0(contrast, "_str"))
  c_grad <- get(paste0("grad_", contrast))

  # set reference to 1 if it is the first level
  reference_n <- which(levels(data[[trt]]) == reference)
  cf_mean_ref <- object$counterfactual.means[reference_n]
  cf_mean_inv <- object$counterfactual.means[3 - reference_n]

  trt_ref <- reference
  trt_inv <- levels(data[[trt]])[-reference_n]

  # apply contrast to point estimate
  marginal_est <- c_est(cf_mean_inv, cf_mean_ref)

  object$marginal_est <- marginal_est
  names(object$marginal_est) <- "marginal_est"
  attr(object$marginal_est, "reference") <- trt_ref
  attr(object$marginal_est, "contrast") <- paste0(contrast, ": ", c_str(trt_inv, trt_ref))

  # apply contrast to varcov
  gr <- c_grad(cf_mean_inv, cf_mean_ref)

  robust_varcov <- object$robust_varcov

  # correct varcov order based on reference
  if (reference_n == 2) diag(robust_varcov) <- rev(diag(robust_varcov))

  marginal_se <- sqrt(t(gr) %*% robust_varcov %*% gr)

  object$marginal_se <- marginal_se[[1]]
  names(object$marginal_se) <- "marginal_se"
  attr(object$marginal_se, "reference") <- trt_ref
  attr(object$marginal_se, "contrast") <- paste0(contrast, ": ", c_str(trt_inv, trt_ref))
  attr(object$marginal_se, "type") <- attr(object$robust_varcov, "type")

  return(object)
}

# risk difference
diff <- function(x, y) {
  x - y
}

diff_str <- function(x, y) {
  paste0(x, "-", y)
}

grad_diff <- function(x, y) {
  c(
    -1,
    1
  )
}

# risk ratio
rr <- function(x, y) {
  x / y
}

rr_str <- function(x, y) {
  paste0(x, "/", y)
}

grad_rr <- function(x, y) {
  c(
    -(x / y^2),
    1 / y
  )
}

# odds ratio
or <- function(x, y) {
  (x / (1 - x)) / (y / (1 - y))
}

or_str <- function(x, y) {
  paste0("odds(", x, ") / odds(", y, ")")
}

grad_or <- function(x, y) {
  # copy of Deriv::Deriv(\(x,y) (x/(1-x)) / (y/(1-y)))
  .e1 <- 1 - y
  .e2 <- 1 - x
  .e3 <- y / .e1
  c(
    -(x * (1 + .e3) / (.e2 * .e1 * .e3^2)),
    .e1 * (1 + x / .e2) / (y * .e2)
  )
}

# log risk ratio
logrr <- function(x, y) {
  log(rr(x, y))
}

logrr_str <- function(x, y) {
  paste0("log(", rr_str(x, y), ")")
}

grad_logrr <- function(x, y) {
  c(
    -(1 / y),
    1 / x
  )
}

# log odds ratio
logor <- function(x, y) {
  log(or(x, y))
}

logor_str <- function(x, y) {
  paste0("log(", or_str(x, y), ")")
}

grad_logor <- function(x, y) {
  c(
    -1 / (y * (1 - y)),
    1 / (x * (1 - x))
  )
}
