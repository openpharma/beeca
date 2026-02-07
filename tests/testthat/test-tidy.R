library(testthat)
library(dplyr)

# Setup test data
trial01$trtp <- factor(trial01$trtp)
trial01 <- trial01 |> dplyr::filter(!is.na(aval))

# Fit basic 2-arm model
fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
  get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

# Fit 3-arm model
data02 <- trial02_cdisc |>
  dplyr::mutate(TRTPN = as.factor(TRTPN))
fit_3arm <- glm(AVAL ~ TRTPN + SEX, family = "binomial", data = data02) |>
  get_marginal_effect(trt = "TRTPN", method = "Ye", contrast = "diff", reference = "1")


# Test basic functionality
test_that("tidy.beeca returns a tibble with expected columns", {
  result <- tidy(fit1)

  expect_s3_class(result, "tbl_df")
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(result)))
})

test_that("tidy.beeca extracts correct number of rows for 2-arm trial", {
  result <- tidy(fit1)

  # Should have 1 row for the single contrast (1 vs 0)
  expect_equal(nrow(result), 1)
})

test_that("tidy.beeca extracts correct number of rows for 3-arm trial", {
  result <- tidy(fit_3arm)

  # Should have 2 rows for two contrasts (default: 2 vs 1, 3 vs 1)
  expect_equal(nrow(result), 2)
})

test_that("tidy.beeca estimates match marginal_est", {
  result <- tidy(fit1)

  # The estimate should match the marginal_est from the beeca object
  expect_equal(result$estimate, as.numeric(fit1$marginal_est), tolerance = 1e-10)
})

test_that("tidy.beeca standard errors match marginal_se", {
  result <- tidy(fit1)

  # The std.error should match the marginal_se from the beeca object
  expect_equal(result$std.error, as.numeric(fit1$marginal_se), tolerance = 1e-10)
})

test_that("tidy.beeca calculates correct test statistics", {
  result <- tidy(fit1)

  # Wald statistic should be estimate / std.error
  expected_stat <- as.numeric(fit1$marginal_est / fit1$marginal_se)
  expect_equal(result$statistic, expected_stat, tolerance = 1e-10)
})

test_that("tidy.beeca calculates correct p-values", {
  result <- tidy(fit1)

  # Two-sided p-value from Wald test
  expected_pval <- 2 * pnorm(-abs(result$statistic))
  expect_equal(result$p.value, expected_pval, tolerance = 1e-10)
})

test_that("tidy.beeca with conf.int adds confidence interval columns", {
  result <- tidy(fit1, conf.int = TRUE)

  expect_true("conf.low" %in% names(result))
  expect_true("conf.high" %in% names(result))
  expect_equal(nrow(result), 1)
})

test_that("tidy.beeca confidence intervals are correctly calculated", {
  result <- tidy(fit1, conf.int = TRUE, conf.level = 0.95)

  z_crit <- qnorm(1 - 0.05 / 2)  # 1.96 for 95% CI
  expected_low <- result$estimate - z_crit * result$std.error
  expected_high <- result$estimate + z_crit * result$std.error

  expect_equal(result$conf.low, expected_low, tolerance = 1e-10)
  expect_equal(result$conf.high, expected_high, tolerance = 1e-10)
})

test_that("tidy.beeca respects conf.level parameter", {
  result_90 <- tidy(fit1, conf.int = TRUE, conf.level = 0.90)
  result_99 <- tidy(fit1, conf.int = TRUE, conf.level = 0.99)

  # 99% CI should be wider than 90% CI
  width_90 <- result_90$conf.high - result_90$conf.low
  width_99 <- result_99$conf.high - result_99$conf.low

  expect_true(width_99 > width_90)
})

test_that("tidy.beeca with include_marginal adds risk estimates", {
  result <- tidy(fit1, include_marginal = TRUE)

  # Should have 3 rows: risk_0, risk_1, and the contrast
  expect_equal(nrow(result), 3)

  # Should have terms starting with "risk_"
  expect_true(any(grepl("^risk_", result$term)))
})

test_that("tidy.beeca marginal risks match counterfactual.means", {
  result <- tidy(fit1, include_marginal = TRUE)

  # Extract risk rows
  risk_rows <- result[grepl("^risk_", result$term), ]

  # Extract treatment levels from term names
  trt_levels <- sub("^risk_", "", risk_rows$term)

  # Check that estimates match counterfactual.means
  for (i in seq_along(trt_levels)) {
    trt_level <- trt_levels[i]
    expect_equal(
      risk_rows$estimate[i],
      as.numeric(fit1$counterfactual.means[trt_level]),
      tolerance = 1e-10
    )
  }
})

test_that("tidy.beeca marginal risk SEs match robust_varcov diagonal", {
  result <- tidy(fit1, include_marginal = TRUE)

  # Extract risk rows
  risk_rows <- result[grepl("^risk_", result$term), ]

  # Extract treatment levels from term names
  trt_levels <- sub("^risk_", "", risk_rows$term)

  # Check that SEs match sqrt of diagonal elements
  for (i in seq_along(trt_levels)) {
    trt_level <- trt_levels[i]
    expect_equal(
      risk_rows$std.error[i],
      sqrt(fit1$robust_varcov[trt_level, trt_level]),
      tolerance = 1e-10
    )
  }
})

test_that("tidy.beeca works with different contrast types", {
  # Test with risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "rr", reference = "0")

  result_rr <- tidy(fit_rr)
  expect_s3_class(result_rr, "tbl_df")
  expect_equal(nrow(result_rr), 1)

  # Test with odds ratio
  fit_or <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "or", reference = "0")

  result_or <- tidy(fit_or)
  expect_s3_class(result_or, "tbl_df")
  expect_equal(nrow(result_or), 1)
})

test_that("tidy.beeca fails gracefully with non-beeca object", {
  regular_glm <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
  skip_if("tidy.glm" %in% methods("tidy"),
          "another package provides tidy.glm, so dispatch succeeds")
  # Without such a method, tidy() has no target for plain glm objects
  expect_error(
    tidy(regular_glm),
    "no applicable method"
  )
})

test_that("tidy.beeca validates conf.int parameter", {
  expect_error(
    tidy(fit1, conf.int = "yes"),
    "conf.int must be logical"
  )

  expect_error(
    tidy(fit1, conf.int = 1),
    "conf.int must be logical"
  )
})

test_that("tidy.beeca validates conf.level parameter", {
  expect_error(
    tidy(fit1, conf.int = TRUE, conf.level = 0),
    "conf.level must be a number between 0 and 1"
  )

  expect_error(
    tidy(fit1, conf.int = TRUE, conf.level = 1),
    "conf.level must be a number between 0 and 1"
  )

  expect_error(
    tidy(fit1, conf.int = TRUE, conf.level = 1.5),
    "conf.level must be a number between 0 and 1"
  )

  expect_error(
    tidy(fit1, conf.int = TRUE, conf.level = -0.5),
    "conf.level must be a number between 0 and 1"
  )
})

test_that("tidy.beeca validates include_marginal parameter", {
  expect_error(
    tidy(fit1, include_marginal = "yes"),
    "include_marginal must be logical"
  )

  expect_error(
    tidy(fit1, include_marginal = 1),
    "include_marginal must be logical"
  )
})

test_that("tidy.beeca handles 3-arm trial with include_marginal", {
  result <- tidy(fit_3arm, include_marginal = TRUE)

  # Should have 3 risk estimates + 2 contrasts = 5 rows
  expect_equal(nrow(result), 5)

  # Should have 3 terms starting with "risk_"
  risk_terms <- result$term[grepl("^risk_", result$term)]
  expect_equal(length(risk_terms), 3)
})

test_that("tidy.beeca p-values are in valid range", {
  result <- tidy(fit1)

  expect_true(all(result$p.value >= 0))
  expect_true(all(result$p.value <= 1))
})

test_that("tidy.beeca term names match TRTVAL in contrasts", {
  result <- tidy(fit1)

  # Extract contrast TRTVAL from marginal_results
  contrast_trtval <- fit1$marginal_results |>
    filter(ANALTYP1 == "INFERENTIAL", !STAT %in% c("risk", "risk_se")) |>
    filter(!grepl("_se$", STAT)) |>
    pull(TRTVAL)

  expect_equal(result$term, contrast_trtval)
})

test_that("tidy.beeca output has no missing values", {
  result <- tidy(fit1, conf.int = TRUE, include_marginal = TRUE)

  expect_false(any(is.na(result)))
})

test_that("tidy.beeca confidence intervals contain the estimate", {
  result <- tidy(fit1, conf.int = TRUE)

  expect_true(all(result$estimate >= result$conf.low))
  expect_true(all(result$estimate <= result$conf.high))
})
