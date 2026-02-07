test_that("format_pvalue formats small p-values correctly", {
  expect_equal(format_pvalue(0.0001), "<0.001")
  expect_equal(format_pvalue(0.05), "0.050")
  expect_equal(format_pvalue(0.001), "0.001")
  expect_equal(format_pvalue(0.0009), "<0.001")
  expect_true(is.na(format_pvalue(NA)))
})

test_that("format_pvalue handles vectors", {
  result <- format_pvalue(c(0.05, 0.001, 0.0001))
  expect_length(result, 3)
  expect_equal(result[1], "0.050")
  expect_equal(result[2], "0.001")
  expect_equal(result[3], "<0.001")
})

test_that("format_pvalue respects custom digits", {
  expect_equal(format_pvalue(0.05, digits = 2), "0.05")
  expect_equal(format_pvalue(0.05, digits = 4), "0.0500")
})

test_that("format_pvalue respects custom threshold", {
  expect_equal(format_pvalue(0.005, small_threshold = 0.01), "<0.01")
  expect_equal(format_pvalue(0.005, small_threshold = 0.001), "0.005")
})

test_that("beeca_summary_table returns expected structure", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  summary_data <- beeca_summary_table(fit)

  expect_type(summary_data, "list")
  expect_named(summary_data, c("arm_statistics", "treatment_effects", "metadata"))

  # Check arm_statistics structure
  expect_s3_class(summary_data$arm_statistics, "tbl_df")
  expect_true("treatment" %in% names(summary_data$arm_statistics))
  expect_true("N" %in% names(summary_data$arm_statistics))
  expect_true("marginal_risk" %in% names(summary_data$arm_statistics))

  # Check treatment_effects structure
  expect_s3_class(summary_data$treatment_effects, "tbl_df")
  expect_true("comparison" %in% names(summary_data$treatment_effects))
  expect_true("estimate" %in% names(summary_data$treatment_effects))
  expect_true("p_value" %in% names(summary_data$treatment_effects))

  # Check metadata
  expect_equal(summary_data$metadata$contrast_type, "diff")
  expect_equal(summary_data$metadata$conf_level, 0.95)
})

test_that("beeca_summary_table handles risk_percent option", {
  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # With percent
  pct <- beeca_summary_table(fit, risk_percent = TRUE)
  # Without percent
  prop <- beeca_summary_table(fit, risk_percent = FALSE)

  # Risk values should differ by factor of 100
  expect_equal(
    pct$arm_statistics$marginal_risk[1],
    prop$arm_statistics$marginal_risk[1] * 100,
    tolerance = 0.001
  )
})

test_that("beeca_summary_table errors on non-beeca object", {
  expect_error(
    beeca_summary_table(lm(mpg ~ wt, data = mtcars)),
    "must be a beeca object"
  )
})

test_that("as_gt errors when gt package not available", {
  skip_if(requireNamespace("gt", quietly = TRUE))

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  expect_error(as_gt(fit), "Package 'gt' is required")
})

test_that("as_gt creates gt table", {
  skip_if_not_installed("gt")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  result <- as_gt(fit)

  expect_s3_class(result, "gt_tbl")
})

test_that("as_gt respects title and subtitle", {
  skip_if_not_installed("gt")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  result <- as_gt(fit, title = "Test Title", subtitle = "Test Subtitle")

  expect_s3_class(result, "gt_tbl")
  # The title should be in the gt object
  expect_true(!is.null(result$`_heading`$title))
})

test_that("as_gt handles different contrasts", {
  skip_if_not_installed("gt")

  trial01$trtp <- factor(trial01$trtp)

  # Risk difference
  fit_diff <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
  expect_s3_class(as_gt(fit_diff), "gt_tbl")

  # Risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
  expect_s3_class(as_gt(fit_rr), "gt_tbl")

  # Odds ratio
  fit_or <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "or", reference = "0")
  expect_s3_class(as_gt(fit_or), "gt_tbl")
})

test_that("as_gt validates inputs", {
  skip_if_not_installed("gt")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Invalid conf.level
  expect_error(as_gt(fit, conf.level = 1.5), "conf.level must be a number between 0 and 1")
  expect_error(as_gt(fit, conf.level = -0.5), "conf.level must be a number between 0 and 1")

  # Non-beeca object
  expect_error(as_gt(lm(mpg ~ wt, data = mtcars)), "must be a beeca object")
})
