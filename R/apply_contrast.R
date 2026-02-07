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
#' @param reference a string or list of strings indicating which treatment
#' group(s) to use as reference level for pairwise comparisons. Accepted values
#' must be a subset of the levels in the treatment variable. Default to the
#' first n-1 treatment levels used in the `glm` object.
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
#' @seealso [predict_counterfactuals()] for generating counterfactual predictions
#' @seealso [average_predictions()] for averaging counterfactual predictions
#' @seealso [estimate_varcov()] for robust variance estimation
#' @seealso [get_marginal_effect()] for estimating marginal effects directly
#' from an original \code{\link[stats]{glm}} object
#' @seealso [beeca_fit()] for streamlined analysis pipeline
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
#' @importFrom utils combn
#' @export
#'
apply_contrast <- function(object, contrast = c("diff", "rr", "or", "logrr", "logor"), reference) {
  # Assert means are available in object
  if (!"counterfactual.means" %in% names(object)) {
    msg <- sprintf(
      "Missing counterfactual means. First run `%1$s <- average_predictions(%1$s)`.",
      deparse(quote(object))
    )
    stop(msg, call. = FALSE)
  }
  # Assert varcov is available in object
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

  # Assert non-missing reference or set default
  if (missing(reference)) {
    reference <- object$xlevels[[trt]][-nlevels(data[[trt]])]
    warning(sprintf("No reference argument was provided, using {%s} as the reference level(s)", paste(reference, collapse = ", ")), call. = FALSE)
  }

  # Assert reference to be subset of the treatment levels
  if (!all(reference %in% levels(data[[trt]]))) {
    stop("Reference levels must be a subset of treatment levels : ", paste(levels(data[[trt]]), collapse = ", "), ".")
  }

  # Only accept at most nlevels(trt)-1 reference levels
  if (length(reference) > (nlevels(data[[trt]])-1)) {
    stop(sprintf("Too many reference levels provided. Expecting at most %s value(s) from the possible treatment levels: {%s}",
                 nlevels(data[[trt]]) - 1,
                 paste(levels(data[[trt]]), collapse = ", ")))
  }

  # Match contrast argument and retrieve relevant helper functions
  contrast <- match.arg(contrast)
  c_est <- get(contrast)
  c_str <- get(paste0(contrast, "_str"))
  c_grad <- get(paste0("grad_", contrast))

  # Extract counterfactual means
  cf_means <- object$counterfactual.means

  # Get all pairwise arm comparisons for contrasts
  # Format according to specified reference levels
  l <- levels(data[[trt]])
  combos <- combn(l, 2, \(x) {
    if (any(reference %in% x)) {
      ref <- intersect(reference, x)[[1]]
      comp <- x[x!=ref]
      c(which(ref == l), which(comp==l))
    } else {
      c(NA, NA)
    }
  })
  combos <- combos[, is.finite(colSums(combos)), drop=F]
  if (length(reference) > 1) {
    combos <- combos[,order(combos[1,])]
  }
  idxs <- matrix(F, nrow=length(l), ncol=length(l))
  rownames(idxs) <- colnames(idxs) <- l
  idxs[cbind(combos[2,], combos[1,])] <- T


  # Get marginal estimates for the contrasts of interest, using correct ref levels
  marginal_est <- outer(cf_means, cf_means, c_est)
  marginal_est <- marginal_est[idxs]

  # Define contrast jacobian
  gr <- matrix(0, nrow=ncol(combos), ncol=length(l))
  # Apply the grad function to all combinations
  contrast_values <- t(apply(combos, 2, function(idx) c_grad(cf_means[idx[2]], cf_means[idx[1]])))
  gr[cbind(1:ncol(combos), combos[1,])] <- contrast_values[, 1]
  gr[cbind(1:ncol(combos), combos[2,])] <- contrast_values[, 2]
  # Get marginal standard error for contrasts of interest
  marginal_se <- sqrt(diag(gr %*% object$robust_varcov %*% t(gr)))

  # Add string description for each comparison
  contrast_desc <- paste0(contrast, ": ",
                          apply(combos, 2, \(x) c_str(l[x[2]], l[x[1]])))
  names(marginal_est) <- contrast_desc
  names(marginal_se) <- contrast_desc

  object$marginal_est <- marginal_est
  attr(object$marginal_est, "reference") <- reference
  attr(object$marginal_est, "contrast") <- contrast_desc

  object$marginal_se <- marginal_se
  attr(object$marginal_se, "reference") <- reference
  attr(object$marginal_se, "contrast") <- contrast_desc
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
