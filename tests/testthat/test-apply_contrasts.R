# set up examples ----------------------------------------------------------

# 2 arm example
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

# 3 arm example
data02 <- trial02_cdisc |>
  transform(TRTPN = as.factor(TRTPN))

fit_3arm <- glm(AVAL ~ TRTPN + SEX, family = "binomial", data = data02) |>
  predict_counterfactuals(trt = "TRTPN") |>
  average_predictions() |>
  estimate_varcov(method="Ye")


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
      apply_contrast(contrast = "diffs", reference = "0"))
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
    "Reference levels must be a subset of treatment levels : A, B.",
    fixed = TRUE
  )

  # Too many reference levels
  expect_error(
    fit |>
      apply_contrast(contrast = "diffs", reference = c("A", "B")),
    "Too many reference levels provided.",
    fixed = TRUE
  )

})


test_that("Check model is sanitized", {
  fit1$sanitized <- NULL
  expect_warning(
    apply_contrast(fit1, reference = "0"),
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


# tests for data with 3 arms ----------------------------------------

test_that("Reference levels correctly handled when >2 treatment levels", {

  expect_warning(
    fit_ref12 <- fit_3arm |> apply_contrast()
  )

  # Default reference levels are first n-1 treatment levels
  expect_equal(
    attr(fit_ref12$marginal_est, "reference"),
    fit_3arm$xlevels$TRTPN[-nlevels(fit_3arm$data$TRTPN)]
  )

  expect_equal(
    attr(fit_ref12$marginal_est, "contrast"),
    c("diff: 2-1", "diff: 3-1", "diff: 3-2")
  )

  # Custom order of reference levels is correctly handled
  fit_ref32 <- fit_3arm |>
    apply_contrast(contrast = "diff", reference = c("3", "2"))

  expect_equal(
    attr(fit_ref32$marginal_est, "contrast"),
    c("diff: 1-2", "diff: 1-3", "diff: 2-3")
  )

  expect_equal(
    fit_ref32$marginal_est[["diff: 1-2"]],
    -fit_ref12$marginal_est[["diff: 2-1"]]
  )
  expect_equal(
    fit_ref32$marginal_est[["diff: 1-3"]],
    -fit_ref12$marginal_est[["diff: 3-1"]]
  )
  expect_equal(
    fit_ref32$marginal_est[["diff: 2-3"]],
    -fit_ref12$marginal_est[["diff: 3-2"]]
  )

  expect_equal(
    unname(fit_ref32$marginal_se)[1:3],
    unname(fit_ref12$marginal_se)[1:3]
  )


  gr <- matrix(c(-1,-1,0,
                 1,0,-1,
                 0,1,1),
               nrow=3, ncol=3)
  expect_equal(
    unname(fit_ref12$marginal_se)[1:3],
    sqrt(diag(gr %*% fit_ref12$robust_varcov %*% t(gr)))
  )


  # Switched references
  fit_ref23 <- fit_3arm |>
    apply_contrast(contrast = "diff", reference = c("2", "3"))

  gr <- matrix(c(1,0,1,
                 -1,-1,0,
                 0,1,-1),
               nrow=3, ncol=3)

  expect_equal(
    unname(fit_ref23$marginal_se)[1:3],
    sqrt(diag(gr %*% fit_ref23$robust_varcov %*% t(gr)))
  )


  # Single reference
  fit_ref2 <- fit_3arm |>
    apply_contrast(contrast = "diff", reference = "2")

  gr <- matrix(c(1,0,-1,
                 -1,0,1),
               nrow=2, ncol=3)
  expect_equal(
    unname(fit_ref2$marginal_se)[1:2],
    sqrt(diag(gr %*% fit_ref2$robust_varcov %*% t(gr)))
  )


  # Check expected values for other contrast type
  # (logor)
  fit_ref12_lor <- fit_3arm |>
    apply_contrast(contrast = "logor", reference = c("1", "2"))

  fit_ref23_lor <- fit_3arm |>
    apply_contrast(contrast = "logor", reference = c("2", "3"))

  expect_equal(
    unname(fit_ref12_lor$marginal_est)[1:3],
    c(1.70521672372, 1.60865867380, -0.09655804992)
  )
  expect_equal(
    unname(fit_ref12_lor$marginal_se)[1:3],
    c(0.3385774965, 0.3336867307, 0.3470551734)
  )

  expect_equal(
    unname(fit_ref23_lor$marginal_est)[1:3],
    c(-1.70521672372, -0.09655804992, -1.60865867380)
  )
  expect_equal(
    unname(fit_ref23_lor$marginal_se)[1:3],
    c(0.3385774965, 0.3470551734, 0.3336867307)
  )

})





# if RobinCar is available compare contrast results
robincar_available <- requireNamespace("RobinCar", quietly = T)

if (robincar_available){

  test_that("Contrast results match RobinCar: diff", {
    rc_diff <- RobinCar::robincar_glm(data.frame(fit_3arm$data),
                                      response_col = "AVAL",
                                      treat_col = "TRTPN",
                                      formula = fit_3arm$formula,
                                      g_family = fit_3arm$family,
                                      contrast_h = "diff")
    fit_ref1 <- fit_3arm |>
      apply_contrast(contrast = "diff", reference = "1")

    expect_equal(
      unname(fit_ref1$marginal_est[1:2]),
      unname(rc_diff$contrast$result$estimate)
    )

    expect_equal(
      unname(fit_ref1$marginal_se[1:2]),
      unname(rc_diff$contrast$result$se)
    )

    # switch reference
    rc_diff_switched <- RobinCar::robincar_glm(data.frame(fit_3arm$data) |>
                                                 transform(TRTPN = factor(TRTPN, levels=c("3", "1", "2"))),
                                               response_col = "AVAL",
                                               treat_col = "TRTPN",
                                               formula = fit_3arm$formula,
                                               g_family = fit_3arm$family,
                                               contrast_h = "diff")
    fit_ref3 <- fit_3arm |>
      apply_contrast(contrast = "diff", reference = "3")

    expect_equal(
      unname(fit_ref3$marginal_est[1:2]),
      unname(rc_diff_switched$contrast$result$estimate)
    )

    expect_equal(
      unname(fit_ref3$marginal_se[1:2]),
      unname(rc_diff_switched$contrast$result$se)
    )

  })


  test_that("Contrast results match RobinCar: risk ratio", {
    rc_rr <- RobinCar::robincar_glm(data.frame(fit_3arm$data),
                                    response_col = "AVAL",
                                    treat_col = "TRTPN",
                                    formula = fit_3arm$formula,
                                    g_family = fit_3arm$family,
                                    contrast_h = "ratio")
    fit_ref1 <- fit_3arm |>
      apply_contrast(contrast = "rr", reference = "1")

    expect_equal(
      unname(fit_ref1$marginal_est[1:2]),
      unname(rc_rr$contrast$result$estimate)
    )

    expect_equal(
      unname(fit_ref1$marginal_se[1:2]),
      unname(rc_rr$contrast$result$se)
    )

    # switch reference
    rc_rr_switched <- RobinCar::robincar_glm(data.frame(fit_3arm$data) |>
                                               transform(TRTPN = factor(TRTPN, levels=c("3", "1", "2"))),
                                             response_col = "AVAL",
                                             treat_col = "TRTPN",
                                             formula = fit_3arm$formula,
                                             g_family = fit_3arm$family,
                                             contrast_h = "ratio")
    fit_ref3 <- fit_3arm |>
      apply_contrast(contrast = "rr", reference = "3")

    expect_equal(
      unname(fit_ref3$marginal_est[1:2]),
      unname(rc_rr_switched$contrast$result$estimate)
    )

    expect_equal(
      unname(fit_ref3$marginal_se[1:2]),
      unname(rc_rr_switched$contrast$result$se)
    )

  })

}
