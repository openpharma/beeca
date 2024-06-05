# set up example ----------------------------------------------------------

data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))

fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
  predict_counterfactuals(trt = "trtp") |>
  average_predictions() |>
  estimate_varcov()

example_varcov <- fit1$robust_varcov
example_theta_s <- fit1$counterfactual.means[["0"]]
example_theta_t <- fit1$counterfactual.means[["1"]]



# test warnings/errors ------------------------------------------------------------


test_that("Assert robust varcov is available in object", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      predict_counterfactuals(trt = "trtp") |>
      average_predictions() |>
      apply_contrast(contrast = "diff", reference = "0"),
    "Missing robust varcov. First run `object <- get_varcov(object, ...)`.",
    fixed = TRUE
  )
})

test_that("Assert counterfactual.means are available in object", {
  fit1[["counterfactual.means"]] <- NULL
  expect_error(
    apply_contrast(object = fit1, contrast = "diff", reference = "0"),
    "Missing counterfactual means. First run `object <- average_predictions(object)`.",
    fixed = TRUE
  )
})

# Confirm calculation of marginal estimates and variances performed
test_that("Test calculation of marginal estimates and variance", {
  data04 <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))

  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data04) |>
    predict_counterfactuals(trt = "trtp") |>
    average_predictions() |>
    estimate_varcov(method = "Ye")

  result <- apply_contrast(fit, contrast = "diff", reference = "1")

  expect_true("marginal_est" %in% names(result))
  expect_true("marginal_se" %in% names(result))
})

# Validate handling of valid and invalid contrast types
test_that("Test valid contrast types", {
  data02 <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data02) |>
    predict_counterfactuals(trt = "trtp") |>
    average_predictions() |>
    estimate_varcov(method = "Ye")
  # Valid contrast
  expect_silent(apply_contrast(fit, contrast = "diff", reference = "1"))
  # Invalid contrast type
  expect_error(apply_contrast(fit, contrast = "unsupported_contrast", reference = "1"))
})

test_that("Correctly throwing errors on mismatched arguments", {
  expect_error(
    fit1 |>
      apply_contrast(contrast = "diffs", reference = "0"),
    "'arg' should be one of"
  )
})

test_that("Test handling of reference levels", {
  data03 <- data.frame(aval = factor(c(1, 0)), trtp = factor(c("A", "B")), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data03) |>
    predict_counterfactuals(trt = "trtp") |>
    average_predictions() |>
    estimate_varcov(method = "Ye")

  # Default reference not explicitly provided, using the first level
  expect_warning(
    apply_contrast(fit, contrast = "diff"),
    "No reference argument was provided"
  )

  # Invalid reference
  expect_error(
    fit |>
      apply_contrast(contrast = "diffs", reference = "C"),
    "Reference must be one of : A, B",
    fixed = TRUE
  )
})


test_that("Check model is sanitized", {
  fit1$sanitized <- NULL
  expect_warning(apply_contrast(fit1, reference = "0"),
    "Input object did not meet the expected format for beeca. Results may not be valid. Consider using ?get_marginal_effect",
    fixed = TRUE
  )
})



# test results ------------------------------------------------------------

## diff ------------------------------------------------------------
test_that("diff contrast is correct for variance estimate", {
  diff_formula <- function(V) {
    V[2, 2] - 2 * V[1, 2] + V[1, 1]
  }

  gr <- grad_diff(example_theta_t, example_theta_s)
  expect_equal(
    (t(gr) %*% example_varcov %*% gr)[[1]],
    diff_formula(example_varcov)
  )

  fit2 <- apply_contrast(fit1, contrast = "diff", reference = "0")
  expect_equal(
    (fit2$marginal_se[[1]])^2,
    diff_formula(example_varcov)
  )
})

test_that("diff contrast is correct for point estimate", {
  expect_equal(
    diff(example_theta_t, example_theta_s),
    example_theta_t - example_theta_s
  )

  fit2 <- apply_contrast(fit1, contrast = "diff", reference = "0")
  expect_equal(
    fit2$marginal_est[[1]],
    example_theta_t - example_theta_s
  )

  # change reference
  fit3 <- apply_contrast(fit1, contrast = "diff", reference = "1")
  expect_equal(
    fit3$marginal_est[[1]],
    example_theta_s - example_theta_t
  )
})


## rr ------------------------------------------------------------
test_that("rr contrast is correct for variance estimate", {
  rr_formula <- function(V, x, y) {
    V[1, 1] * (x^2 / y^4) - (V[1, 2] * x) / (y^3) - (V[2, 1] * x) / (y^3) + V[2, 2] / (y^2)
  }

  gr <- grad_rr(example_theta_t, example_theta_s)
  expect_equal(
    (t(gr) %*% example_varcov %*% gr)[[1]],
    rr_formula(example_varcov, example_theta_t, example_theta_s)
  )

  fit2 <- apply_contrast(fit1, contrast = "rr", reference = "0")
  expect_equal(
    (fit2$marginal_se[[1]])^2,
    rr_formula(example_varcov, example_theta_t, example_theta_s)
  )

  # change reference
  fit3 <- apply_contrast(fit1, contrast = "rr", reference = "1")
  expect_equal(
    (fit3$marginal_se[[1]])^2,
    rr_formula(
      matrix(rev(example_varcov), 2, 2),
      example_theta_s, example_theta_t
    )
  )
})


## logrr ------------------------------------------------------------
test_that("logrr contrast is correct for variance estimate", {
  logrr_formula <- function(V, x, y) {
    V[2, 2] / (x^2) - (2 * V[1, 2]) / (x * y) + V[1, 1] / (y^2)
  }

  gr <- grad_logrr(example_theta_t, example_theta_s)
  expect_equal(
    (t(gr) %*% example_varcov %*% gr)[[1]],
    logrr_formula(example_varcov, example_theta_t, example_theta_s)
  )

  fit2 <- apply_contrast(fit1, contrast = "logrr", reference = "0")
  expect_equal(
    (fit2$marginal_se[[1]])^2,
    logrr_formula(example_varcov, example_theta_t, example_theta_s)
  )

  # change reference
  fit3 <- apply_contrast(fit1, contrast = "logrr", reference = "1")
  expect_equal(
    (fit3$marginal_se[[1]])^2,
    logrr_formula(
      matrix(rev(example_varcov), 2, 2),
      example_theta_s, example_theta_t
    )
  )
})


## or ------------------------------------------------------------
test_that("or contrast is correct for variance estimate", {
  or_formula <- function(V, x, y) {
    gx <- -(x * (1 + (y / (1 - y)))) / ((1 - x) * (1 - y) * (y / (1 - y))^2)
    gy <- (1 - y) * (1 + (x / (1 - x))) / (y * (1 - x))

    V[1, 1] * (gx^2) + V[1, 2] * gx * gy + V[2, 1] * gx * gy + V[2, 2] * (gy^2)
  }

  gr <- grad_or(example_theta_t, example_theta_s)
  expect_equal(
    (t(gr) %*% example_varcov %*% gr)[[1]],
    or_formula(example_varcov, example_theta_t, example_theta_s)
  )

  fit2 <- apply_contrast(fit1, contrast = "or", reference = "0")
  expect_equal(
    (fit2$marginal_se[[1]])^2,
    or_formula(example_varcov, example_theta_t, example_theta_s)
  )


  # change reference
  fit3 <- apply_contrast(fit1, contrast = "or", reference = "1")
  expect_equal(
    (fit3$marginal_se[[1]])^2,
    or_formula(
      matrix(rev(example_varcov), 2, 2),
      example_theta_s, example_theta_t
    )
  )
})


## logor ------------------------------------------------------------
test_that("logor contrast is correct for variance estimate", {
  logor_formula <- function(V, x, y) {
    V[2, 2] / ((x^2) * (1 - x)^2) - (2 * V[1, 2]) / (x * (1 - x) * y * (1 - y)) + V[1, 1] / ((y^2) * (1 - y)^2)
  }

  gr <- grad_logor(example_theta_t, example_theta_s)
  expect_equal(
    (t(gr) %*% example_varcov %*% gr)[[1]],
    logor_formula(example_varcov, example_theta_t, example_theta_s)
  )

  fit2 <- apply_contrast(fit1, contrast = "logor", reference = "0")
  expect_equal(
    (fit2$marginal_se[[1]])^2,
    logor_formula(example_varcov, example_theta_t, example_theta_s)
  )


  # change reference
  fit3 <- apply_contrast(fit1, contrast = "logor", reference = "1")
  expect_equal(
    (fit3$marginal_se[[1]])^2,
    logor_formula(
      matrix(rev(example_varcov), 2, 2),
      example_theta_s, example_theta_t
    )
  )
})
