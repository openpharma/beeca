test_that("print.beeca displays basic output correctly", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Capture output
  output <- capture.output(print(fit))
  output_text <- paste(output, collapse = "\n")

  # Check key components are present
  expect_true(grepl("beeca:", output_text))
  expect_true(grepl("Contrast:", output_text))
  expect_true(grepl("Estimate:", output_text))
  expect_true(grepl("SE =", output_text))
  expect_true(grepl("Z-value:", output_text))
  expect_true(grepl("P-value:", output_text))
})

test_that("print.beeca handles multiple contrasts", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff")

  output <- capture.output(print(fit))
  output_text <- paste(output, collapse = "\n")

  # Should show multiple comparisons
  expect_true(grepl("Contrast:", output_text))
})

test_that("print.beeca respects digits parameter", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Test with different digits
  output_2 <- capture.output(print(fit, digits = 2))
  output_6 <- capture.output(print(fit, digits = 6))

  # Outputs should differ in precision
  expect_false(identical(output_2, output_6))
})

test_that("print.beeca shows helpful usage hints", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  output <- capture.output(print(fit))
  output_text <- paste(output, collapse = "\n")

  # Check for helpful hints
  expect_true(grepl("summary\\(\\)", output_text) ||
                grepl("tidy\\(\\)", output_text) ||
                grepl("plot\\(\\)", output_text))
})

test_that("print.beeca returns object invisibly", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # print should return the object invisibly
  result <- withVisible(print(fit))
  expect_false(result$visible)
  expect_identical(result$value, fit)
})

test_that("summary.beeca displays comprehensive output", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  # Check main sections are present
  expect_true(grepl("Covariate-Adjusted Marginal Treatment Effect Analysis", output_text))
  expect_true(grepl("Model Information:", output_text))
  expect_true(grepl("Marginal Risks", output_text))
  expect_true(grepl("Treatment Effect Estimates:", output_text))
})

test_that("summary.beeca shows model information correctly", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  # Check for model details
  expect_true(grepl("Working model:", output_text))
  expect_true(grepl("Variance estimator:", output_text))
  expect_true(grepl("Ye", output_text))
  expect_true(grepl("Contrast type:", output_text))
  expect_true(grepl("Sample size:", output_text))
  expect_true(grepl("Outcome variable:", output_text))
})

test_that("summary.beeca shows marginal risks with CI", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  # Check for risk table components
  expect_true(grepl("Treatment", output_text))
  expect_true(grepl("Risk", output_text))
  expect_true(grepl("SE", output_text))
  expect_true(grepl("95% CI", output_text))
})

test_that("summary.beeca shows treatment effects with statistics", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  # Check for effects table components
  expect_true(grepl("Comparison", output_text))
  expect_true(grepl("Estimate", output_text))
  expect_true(grepl("Z value", output_text))
  expect_true(grepl("P value", output_text))
})

test_that("summary.beeca respects conf.level parameter", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Test with 90% CI
  output_90 <- capture.output(summary(fit, conf.level = 0.90))
  output_text_90 <- paste(output_90, collapse = "\n")
  expect_true(grepl("90% CI", output_text_90))

  # Test with 99% CI
  output_99 <- capture.output(summary(fit, conf.level = 0.99))
  output_text_99 <- paste(output_99, collapse = "\n")
  expect_true(grepl("99% CI", output_text_99))
})

test_that("summary.beeca respects digits parameter", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Different digits should produce different output
  output_2 <- capture.output(summary(fit, digits = 2))
  output_6 <- capture.output(summary(fit, digits = 6))

  expect_false(identical(output_2, output_6))
})

test_that("summary.beeca handles multiple comparisons", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  # Should show all comparisons
  expect_true(grepl("Comparison", output_text))
})

test_that("summary.beeca shows correct contrast type labels", {
  trial01$trtp <- factor(trial01$trtp)

  # Risk difference
  fit_diff <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
  output_diff <- capture.output(summary(fit_diff))
  expect_true(any(grepl("diff", output_diff)))

  # Risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
  output_rr <- capture.output(summary(fit_rr))
  expect_true(any(grepl("rr", output_rr)))
})

test_that("summary.beeca returns object invisibly", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  result <- withVisible(summary(fit))
  expect_false(result$visible)
  expect_identical(result$value, fit)
})

test_that("summary.beeca handles Ge method correctly", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "diff", reference = "0")

  output <- capture.output(summary(fit))
  output_text <- paste(output, collapse = "\n")

  expect_true(grepl("Ge", output_text))
})

test_that("print.beeca handles different contrast types", {
  trial01$trtp <- factor(trial01$trtp)

  # Test risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
  output_rr <- capture.output(print(fit_rr))
  expect_true(any(grepl("rr:", output_rr)))

  # Test odds ratio
  fit_or <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "or", reference = "0")
  output_or <- capture.output(print(fit_or))
  expect_true(any(grepl("or:", output_or)))
})
