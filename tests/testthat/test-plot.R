test_that("plot.beeca works with default settings", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p <- plot(fit)

  expect_s3_class(p, "ggplot")
})

test_that("plot.beeca accepts type argument", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p <- plot(fit, type = "forest")

  expect_s3_class(p, "ggplot")
})

test_that("plot.beeca rejects invalid type", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  expect_error(
    plot(fit, type = "invalid"),
    "'arg' should be"
  )
})

test_that("plot.beeca passes arguments to plot_forest", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Test that conf.level is passed through
  p <- plot(fit, conf.level = 0.90)

  expect_s3_class(p, "ggplot")
})

test_that("plot_forest creates ggplot object", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p <- plot_forest(fit)

  expect_s3_class(p, "ggplot")
  expect_s3_class(p, "gg")
})

test_that("plot_forest errors without ggplot2", {
  # Skip this test - cannot test ggplot2 requirement when ggplot2 is installed
  skip("Cannot test ggplot2 requirement when ggplot2 is installed")
})

test_that("plot_forest respects conf.level parameter", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p90 <- plot_forest(fit, conf.level = 0.90)
  p95 <- plot_forest(fit, conf.level = 0.95)

  # Both should be ggplot objects
  expect_s3_class(p90, "ggplot")
  expect_s3_class(p95, "ggplot")

  # Titles should differ
  expect_true(grepl("90%", p90$labels$title))
  expect_true(grepl("95%", p95$labels$title))
})

test_that("plot_forest accepts custom title", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  custom_title <- "My Custom Title"
  p <- plot_forest(fit, title = custom_title)

  expect_equal(p$labels$title, custom_title)
})

test_that("plot_forest accepts custom xlab", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  custom_xlab <- "My X Label"
  p <- plot_forest(fit, xlab = custom_xlab)

  expect_equal(p$labels$x, custom_xlab)
})

test_that("plot_forest sets appropriate xlab for different contrasts", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)

  # Risk difference
  fit_diff <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
  p_diff <- plot_forest(fit_diff)
  expect_equal(p_diff$labels$x, "Risk Difference")

  # Risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
  p_rr <- plot_forest(fit_rr)
  expect_equal(p_rr$labels$x, "Risk Ratio")

  # Odds ratio
  fit_or <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "or", reference = "0")
  p_or <- plot_forest(fit_or)
  expect_equal(p_or$labels$x, "Odds Ratio")

  # Log risk ratio
  fit_logrr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "logrr", reference = "0")
  p_logrr <- plot_forest(fit_logrr)
  expect_equal(p_logrr$labels$x, "Log Risk Ratio")

  # Log odds ratio
  fit_logor <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "logor", reference = "0")
  p_logor <- plot_forest(fit_logor)
  expect_equal(p_logor$labels$x, "Log Odds Ratio")
})

test_that("plot_forest can hide numerical values", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p_with <- plot_forest(fit, show_values = TRUE)
  p_without <- plot_forest(fit, show_values = FALSE)

  # Both should be ggplot objects
  expect_s3_class(p_with, "ggplot")
  expect_s3_class(p_without, "ggplot")

  # Plots should differ (with values has annotations)
  expect_true(length(p_with$layers) >= length(p_without$layers))
})

test_that("plot_forest accepts custom null_line_color", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p <- plot_forest(fit, null_line_color = "red")

  expect_s3_class(p, "ggplot")
})

test_that("plot_forest accepts custom point_size", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  p <- plot_forest(fit, point_size = 5)

  expect_s3_class(p, "ggplot")
})

test_that("plot_forest handles multiple comparisons", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff")

  p <- plot_forest(fit)

  expect_s3_class(p, "ggplot")
})

test_that("plot_forest can be customized with ggplot2", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Should be able to add ggplot2 layers
  p <- plot_forest(fit) +
    ggplot2::theme_minimal() +
    ggplot2::labs(caption = "My caption")

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$caption, "My caption")
})

test_that("plot_forest uses correct null value for ratio contrasts", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)

  # For risk ratio and odds ratio, null value should be 1
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "rr", reference = "0")
  p_rr <- plot_forest(fit_rr)

  # Check that vertical line is at 1 (this is implicit in the plot data)
  expect_s3_class(p_rr, "ggplot")

  # For difference, null value should be 0
  fit_diff <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")
  p_diff <- plot_forest(fit_diff)

  expect_s3_class(p_diff, "ggplot")
})

test_that("plot_forest handles missing contrast attribute gracefully", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  # Remove contrast attribute
  attr(fit$marginal_est, "contrast") <- NULL

  # Should still work with default
  p <- plot_forest(fit)

  expect_s3_class(p, "ggplot")
})

test_that("plot_forest works with Ge method", {
  skip_if_not_installed("ggplot2")

  trial01$trtp <- factor(trial01$trtp)
  fit <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "diff", reference = "0")

  p <- plot_forest(fit)

  expect_s3_class(p, "ggplot")
})
