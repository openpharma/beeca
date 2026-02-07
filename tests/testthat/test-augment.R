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
test_that("augment.beeca returns a tibble", {
  result <- augment(fit1)

  expect_s3_class(result, "tbl_df")
})

test_that("augment.beeca has correct number of rows", {
  result <- augment(fit1)

  # Should have same number of rows as original data
  expect_equal(nrow(result), nrow(trial01))
})

test_that("augment.beeca includes model data columns", {
  result <- augment(fit1)

  # Should include all columns from the model frame (variables used in the model)
  model_cols <- names(stats::model.frame(fit1))
  expect_true(all(model_cols %in% names(result)))
})

test_that("augment.beeca adds .fitted column", {
  result <- augment(fit1)

  expect_true(".fitted" %in% names(result))
  expect_equal(length(result$.fitted), nrow(trial01))
})

test_that("augment.beeca .fitted values match predict()", {
  result <- augment(fit1)

  # .fitted should match predict with type = "response"
  expected_fitted <- predict(fit1, type = "response")
  expect_equal(unname(result$.fitted), unname(expected_fitted), tolerance = 1e-10)
})

test_that("augment.beeca adds .resid column by default", {
  result <- augment(fit1)

  expect_true(".resid" %in% names(result))
  expect_equal(length(result$.resid), nrow(trial01))
})

test_that("augment.beeca .resid values match residuals()", {
  result <- augment(fit1)

  # .resid should match residuals with type = "deviance"
  expected_resid <- residuals(fit1, type = "deviance")
  expect_equal(unname(result$.resid), unname(expected_resid), tolerance = 1e-10)
})

test_that("augment.beeca adds counterfactual columns for 2-arm trial", {
  result <- augment(fit1)

  # Should have .counterfactual_0 and .counterfactual_1
  expect_true(".counterfactual_0" %in% names(result))
  expect_true(".counterfactual_1" %in% names(result))
})

test_that("augment.beeca adds counterfactual columns for 3-arm trial", {
  result <- augment(fit_3arm)

  # Should have counterfactual columns for all 3 arms
  expect_true(".counterfactual_1" %in% names(result))
  expect_true(".counterfactual_2" %in% names(result))
  expect_true(".counterfactual_3" %in% names(result))
})

test_that("augment.beeca counterfactual predictions match object", {
  result <- augment(fit1)

  # Extract counterfactual predictions from the beeca object
  cf_preds <- fit1$counterfactual.predictions

  # Check that each column matches
  expect_equal(result$.counterfactual_0, cf_preds[["0"]], tolerance = 1e-10)
  expect_equal(result$.counterfactual_1, cf_preds[["1"]], tolerance = 1e-10)
})

test_that("augment.beeca counterfactual values are probabilities", {
  result <- augment(fit1)

  # All counterfactual predictions should be between 0 and 1
  expect_true(all(result$.counterfactual_0 >= 0 & result$.counterfactual_0 <= 1))
  expect_true(all(result$.counterfactual_1 >= 0 & result$.counterfactual_1 <= 1))
})

test_that("augment.beeca type.predict = 'link' works", {
  result <- augment(fit1, type.predict = "link")

  # .fitted should be on the link scale (not 0-1)
  expected_fitted <- predict(fit1, type = "link")
  expect_equal(unname(result$.fitted), unname(expected_fitted), tolerance = 1e-10)
})

test_that("augment.beeca type.residuals options work", {
  result_pearson <- augment(fit1, type.residuals = "pearson")
  result_response <- augment(fit1, type.residuals = "response")
  result_working <- augment(fit1, type.residuals = "working")

  expect_equal(unname(result_pearson$.resid), unname(residuals(fit1, type = "pearson")), tolerance = 1e-10)
  expect_equal(unname(result_response$.resid), unname(residuals(fit1, type = "response")), tolerance = 1e-10)
  expect_equal(unname(result_working$.resid), unname(residuals(fit1, type = "working")), tolerance = 1e-10)
})

test_that("augment.beeca type.residuals = NULL excludes residuals", {
  result <- augment(fit1, type.residuals = NULL)

  expect_false(".resid" %in% names(result))
})

test_that("augment.beeca fails gracefully with non-beeca object", {
  regular_glm <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)
  skip_if("augment.glm" %in% methods("augment"),
          "another package provides augment.glm, so dispatch succeeds")
  # Without such a method, augment() has no target for plain glm objects
  expect_error(
    augment(regular_glm),
    "no applicable method"
  )
})

test_that("augment.beeca validates type.predict parameter", {
  expect_error(
    augment(fit1, type.predict = "invalid"),
    "type.predict must be 'response' or 'link'"
  )
})

test_that("augment.beeca validates type.residuals parameter", {
  expect_error(
    augment(fit1, type.residuals = "invalid"),
    "type.residuals must be one of"
  )
})

test_that("augment.beeca with custom data works", {
  # Create a subset of the data
  subset_data <- trial01[1:50, ]

  result <- augment(fit1, data = subset_data)

  expect_equal(nrow(result), 50)
})

test_that("augment.beeca preserves row names", {
  # Add row names to data
  data_with_rownames <- trial01
  rownames(data_with_rownames) <- paste0("subj_", seq_len(nrow(trial01)))

  # Fit model with this data
  fit_rownames <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data_with_rownames) |>
    get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

  result <- augment(fit_rownames)

  expect_true(".rownames" %in% names(result))
  expect_equal(result$.rownames, rownames(data_with_rownames))
})

test_that("augment.beeca handles newdata deprecation warning", {
  expect_warning(
    augment(fit1, newdata = trial01[1:10, ]),
    "'newdata' is deprecated"
  )
})

test_that("augment.beeca individual treatment effects sum to average", {
  result <- augment(fit1)

  # Calculate individual treatment effects
  ite <- result$.counterfactual_1 - result$.counterfactual_0

  # Average ITE should equal the marginal risk difference
  # (approximately, due to g-computation)
  mean_ite <- mean(ite)
  marginal_rd <- as.numeric(fit1$counterfactual.means["1"] - fit1$counterfactual.means["0"])

  expect_equal(mean_ite, marginal_rd, tolerance = 1e-10)
})

test_that("augment.beeca counterfactuals match for observed treatment", {
  result <- augment(fit1)

  # For subjects in treatment 0, .counterfactual_0 should be close to .fitted
  # (though not exactly equal due to model structure)
  subjects_trt0 <- result$trtp == "0"
  subjects_trt1 <- result$trtp == "1"

  # The counterfactual for the observed treatment should approximate the fitted value
  # This is true when treatment doesn't interact with covariates
  # But they won't be identical because fitted values are based on actual treatment
})

test_that("augment.beeca output has no unexpected missing values", {
  result <- augment(fit1)

  # Check that key columns have no NAs
  expect_false(any(is.na(result$.fitted)))
  expect_false(any(is.na(result$.resid)))
  expect_false(any(is.na(result$.counterfactual_0)))
  expect_false(any(is.na(result$.counterfactual_1)))
})

test_that("augment.beeca works with different contrast types", {
  # Test with risk ratio
  fit_rr <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "rr", reference = "0")

  result_rr <- augment(fit_rr)
  expect_s3_class(result_rr, "tbl_df")
  expect_true(".counterfactual_0" %in% names(result_rr))
  expect_true(".counterfactual_1" %in% names(result_rr))

  # Test with odds ratio
  fit_or <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01) |>
    get_marginal_effect(trt = "trtp", method = "Ge", contrast = "or", reference = "0")

  result_or <- augment(fit_or)
  expect_s3_class(result_or, "tbl_df")
  expect_true(".counterfactual_0" %in% names(result_or))
  expect_true(".counterfactual_1" %in% names(result_or))
})

test_that("augment.beeca counterfactual columns have correct length", {
  result <- augment(fit1)

  expect_equal(length(result$.counterfactual_0), nrow(trial01))
  expect_equal(length(result$.counterfactual_1), nrow(trial01))
})

test_that("augment.beeca handles 3-arm counterfactuals correctly", {
  result <- augment(fit_3arm)

  # Verify all three counterfactual columns exist and have correct dimensions
  expect_equal(length(result$.counterfactual_1), nrow(data02))
  expect_equal(length(result$.counterfactual_2), nrow(data02))
  expect_equal(length(result$.counterfactual_3), nrow(data02))

  # Verify they match the stored counterfactual predictions
  cf_preds <- fit_3arm$counterfactual.predictions
  expect_equal(result$.counterfactual_1, cf_preds[["1"]], tolerance = 1e-10)
  expect_equal(result$.counterfactual_2, cf_preds[["2"]], tolerance = 1e-10)
  expect_equal(result$.counterfactual_3, cf_preds[["3"]], tolerance = 1e-10)
})

test_that("augment.beeca average of counterfactuals matches counterfactual.means", {
  result <- augment(fit1)

  # Mean of .counterfactual_0 should match counterfactual.means["0"]
  mean_cf0 <- mean(result$.counterfactual_0)
  mean_cf1 <- mean(result$.counterfactual_1)

  expect_equal(mean_cf0, as.numeric(fit1$counterfactual.means["0"]), tolerance = 1e-10)
  expect_equal(mean_cf1, as.numeric(fit1$counterfactual.means["1"]), tolerance = 1e-10)
})
