#' (internal) Sanitize functions to check model and data within GLM model object
#'
#' @description
#'
#' Performs checks on a GLM model object to ensure it meets specific criteria
#' required for further analysis using other functions from the `beeca` package.
#'
#' This includes verifying the model's family, link function, data completeness
#' and mode convergence.
#'
#' Currently it supports models with a binomial family and canonical logit link.
#'
#' @param model a model object, currently only
#' \link[stats]{glm} with binomial family canonical link is supported.
#' @param ... arguments passed to or from other methods.
#' @return if model is non-compliant will throw warnings or errors.
#' @keywords internal
#' @export
sanitize_model <- function(model, ...) {
  UseMethod("sanitize_model")
}

#' (internal) Sanitize a glm model
#' @param model a \link[stats]{glm} with binomial family canonical link.
#' @param trt the name of the treatment variable on the right-hand side of the formula in a \link[stats]{glm}.
#' @param ... ignored.
#' @return if model is non-compliant will throw warnings or errors.
#' @importFrom stats model.frame model.matrix terms
#' @keywords internal
#' @export
#' @examples
#' trial01$trtp <- factor(trial01$trtp)
#' fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
#' fit1 <- sanitize_model(fit1, "trtp")
#'
sanitize_model.glm <- function(model, trt, ...) {
  # sanitize variable
  sanitize_variable(model, trt)

  # sanitize model

  reasons_stop <- reasons_warn <- NULL

  # check family and link function
  if (model$family$family != "binomial" | model$family$link != "logit") {
    reasons_stop <- c(reasons_stop, "not in the binomial family with logit link function")
  }

  # check interaction
  if (any(attr(model$terms, "order") > 1)) {
    interactions <- attr(model$terms, "term.labels")[which(attr(model$terms, "order") > 1)]
    if (trt %in% unlist(strsplit(interactions, ":"))) {
      reasons_stop <- c(reasons_stop, "with treatment-covariate interaction terms")
    }
  }

  # check for missing data
  if (nrow(model$data) != nrow(stats::model.frame(model))) {
    n_missing <- nrow(model$data) - nrow(model.frame(model))
    if (n_missing == 1) {
      reasons_warn <- c(reasons_warn, "There is 1 record omitted from the original data due to missing values, please check if they should be imputed prior to model fitting.")
    } else {
      reasons_warn <- c(reasons_warn, sprintf("There are %s records omitted from the original data due to missing values, please check if they should be imputed prior to model fitting.", n_missing))
    }
  }

  # check for full rank
  if (ncol(stats::model.matrix(model)) != model$qr$rank) {
    reasons_stop <- c(reasons_stop, "The data does not have full rank, please check glm model fitting.")
  }

  # check convergence
  if (!model$converged) {
    reasons_warn <- c(reasons_warn, "The glm model was not converged, please check glm model fitting.")
  }

  # print error messages
  if (!is.null(reasons_stop)) {
    msg_stop <- sapply(reasons_stop, function(x) c(sprintf("Model of class glm %s is not supported.", x), "\n"))
    stop(msg_stop, call. = FALSE)
  }

  # print warning messages
  if (!is.null(reasons_warn)) {
    for (msg_warn in reasons_warn) warning(msg_warn, call. = FALSE)
  }

  model$sanitized <- TRUE

  return(model)
}

#' @export
sanitize_model.default <- function(model, trt, ...) {
  if (!inherits(model, "glm")) {
    msg <- c(sprintf('Model of class "%s" is not supported.', class(model)[1]))
    stop(msg, call. = FALSE)
  }
}

#' (internal) Sanitize function to check model and data
#' @param model an \link[stats]{glm} model object.
#' @param trt the name of the treatment variable on the right-hand side of the glm formula.
#' @return if model and variable are non-compliant, will throw warnings or error.
#' @keywords internal
sanitize_variable <- function(model, trt) {
  data <- .get_data(model)

  # check trt is on right hand side of model formula
  if (!trt %in% attr(model$terms, "term.labels")) {
    msg <- sprintf('Did not find the treatment variable "%s" on right hand side of the model formula "%s".', trt, format(model$formula))
    stop(msg, call. = FALSE)
  }

  # assert trt is a factor
  if (!is.factor(data[[trt]])) {
    msg <- sprintf('Treatment variable "%s" must be of type factor, not "%s".', trt, typeof(data[[trt]]))
    stop(msg, call. = FALSE)
  }
  # assert trt has 2 levels
  if (nlevels(data[[trt]]) != 2) {
    msg <- sprintf('Treatment variable "%s" must have 2 levels. Found %s: {%s}.', trt, length(levels(data[[trt]])), paste(levels(data[[trt]]), collapse = ","))
    stop(msg, call. = FALSE)
  }
  # assert response is coded 0/1
  response_var <- with(attributes(terms(model)), as.character(variables[response + 1]))
  if (!identical(levels(factor(data[[response_var]])), c("0", "1"))) {
    msg <- sprintf('Response variable must be coded as "0" (no event) / "1" (event). Found %s.', paste(levels(factor(data[[response_var]])), collapse = " / "))
    stop(msg, call. = FALSE)
  }
  # # Inform user of internal treatment level coding
  # if (levels(data[[trt]])[1] != 0 | levels(data[[trt]])[2] != 1){
  #   message(sprintf("Treatment level coding:\n\t%s = 0\n\t%s = 1", levels(data[[trt]])[1], levels(data[[trt]])[2]))
  # }
}

.get_data <- function(model) {
  return(model$model)
}

.assert_sanitized <- function(object, trt, warn = F) {
  if (is.null(object$sanitized) || object$sanitized == FALSE) {
    if (warn) {
      warning("Input object did not meet the expected format for beeca. Results may not be valid. Consider using ?get_marginal_effect",
        call. = FALSE
      )
    }
    object <- sanitize_model(object, trt)
  }

  return(object)
}
