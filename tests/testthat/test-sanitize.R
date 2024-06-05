# set up example ----------------------------------------------------------

data <- trial01
data$trtp <- factor(data$trtp)
data_complete <- na.omit(data)

# Correct model type and properties
test_that("Correct model passes with no errors or warnings", {
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data_complete)
  expect_silent(sanitize_model(fit, "trtp"))
})

# test warnings/errors ------------------------------------------------------------

test_that("Correctly throwing warnings on missing value", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data)
  expect_warning(
    sanitize_model(fit1, "trtp"),
    "There is 1 record omitted from the original data due to missing values, please check if they should be imputed prior to model fitting."
  )

  data$aval[10] <- NA
  fit2 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data)
  expect_warning(
    sanitize_model(fit2, "trtp"),
    "There are 2 records omitted from the original data due to missing values, please check if they should be imputed prior to model fitting."
  )
})

test_that("Correctly throwing errors on incompatible link function", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = binomial(link = "probit"), data = data)
  expect_error(
    sanitize_model(fit1, "trtp"),
    "Model of class glm not in the binomial family with logit link function is not supported."
  )
})

test_that("Correctly throwing errors on treatment interaction", {
  fit1 <- glm(aval ~ trtp * bl_cov, family = "binomial", data = data)
  expect_error(
    sanitize_model(fit1, "trtp"),
    "Model of class glm with treatment-covariate interaction terms is not supported."
  )
})

# Treatment variable in model data
test_that("Treatment variable in model data", {
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data_complete)
  expect_error(
    sanitize_model(fit, "trt"),
    'Did not find the treatment variable "trt" on right hand side of the model formula',
  )
})

test_that("Check treatment variable is a factor", {
  data_complete$trtp <- as.numeric(data_complete$trtp)
  fit1 <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data_complete)
  expect_error(sanitize_model(fit1, "trtp"),
    'Treatment variable "trtp" must be of type factor, not "double".',
    fixed = TRUE
  )
})

test_that("Correctly throwing errors on treatment variable > 2 levels", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01 |> transform(trtp = factor(as.numeric(trtp) + rbinom(268, size = 1, prob = 0.5))))
  expect_error(sanitize_model(fit1, "trtp"),
    'Treatment variable "trtp" must have 2 levels. Found 3: {0,1,2}.',
    fixed = TRUE
  )
})

test_that("Check response variable is 0/1", {
  data_complete$aval <- replace(data_complete$aval, data_complete$aval == "0", "0.5")
  levels(data_complete[["aval"]]) <- c("0.5", "1")
  data_complete$aval <- as.factor(data_complete$aval)
  fit1 <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data_complete)
  expect_error(sanitize_model(fit1, "trtp"))
})

test_that("Correctly throw error on incorrect model class", {
  lm_fit <- stats::lm(aval ~ trtp + bl_cov, data = data_complete)

  expect_error(sanitize_model(lm_fit, "trtp"),
    'Model of class "lm" is not supported.',
    fixed = TRUE
  )
})


test_that("Throw warning if model matrix not full rank", {
  # create rank deficient example
  mat <- data.frame(list(
    y = rbinom(100, 1, 0.5),
    trtp = factor(rbinom(100, 1, 0.5)),
    x1 = rnorm(100),
    x2 = rnorm(100)
  ))
  mat[["x3"]] <- mat$x1 + mat$x2
  fit1 <- glm(y ~ trtp + x1 + x2 + x3, family = "binomial", data = mat)

  expect_error(sanitize_model(fit1, "trtp"),
    "The data does not have full rank, please check glm model fitting.",
    fixed = TRUE
  )
})


test_that("Throw warning if model not converged", {
  # fit glm with reduced max iterations so does not converge
  suppressWarnings(
    fit1 <- glm(aval ~ trtp,
      family = "binomial", data = data,
      control = glm.control(maxit = 1)
    )
  )

  expect_warning(sanitize_model(fit1, "trtp"),
    "The glm model was not converged, please check glm model fitting.",
    fixed = TRUE
  )
})


test_that("Throw error if treatment not on right hand side of model formula", {
  fit1 <- glm(trtp ~ aval, family = "binomial", data = data)

  expect_error(sanitize_variable(fit1, "trtp"),
    'Did not find the treatment variable "trtp" on right hand side of the model formula',
    fixed = TRUE
  )
})
