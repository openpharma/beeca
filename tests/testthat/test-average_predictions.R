# set up example ----------------------------------------------------------

data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))

# test calculation --------------------------------------------------------

test_that("averaging counterfactual predictions works", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01)
  fit1 <- predict_counterfactuals(fit1, "trtp")
  fit1 <- average_predictions(fit1)

  expected_output <- c(
    "0" = mean(fit1$counterfactual.predictions$cf_pred_0),
    "1" = mean(fit1$counterfactual.predictions$cf_pred_1)
  )

  expect_equal(fit1$counterfactual.means, expected_output)
})

test_that("averaging throws error if no counterfactual predictions in glm object", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01)
  expect_error(average_predictions(fit1))
})



# Test Case 1: Missing Counterfactual Predictions
test_that("Test missing counterfactual predictions", {
  data <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data)
  expect_error(
    average_predictions(fit),
    "Missing counterfactual predictions"
  )
})

# Test Case 2: Correct Calculation of Averages
test_that("Test correct calculation of averages", {
  data <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data)
  fit1 <- predict_counterfactuals(fit, "trtp")
  result <- average_predictions(fit1)
  computed_means <- colMeans(fit1$counterfactual.predictions)
  names(computed_means) <- levels(data$trtp)
  expect_equal(result$counterfactual.means, computed_means)
})

# Test Case 3: Dependency on Predict Counterfactuals
test_that("Ensure function depends on predict_counterfactuals", {
  data <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data)
  fit$counterfactual.predictions <- matrix(runif(200), nrow = 100, ncol = 2)
  # Assuming predict_counterfactuals modifies object in a specific way that average_predictions expects
  fit <- predict_counterfactuals(fit, trt = "trtp")
  expect_silent(average_predictions(fit))
})

# Test Case 4: Integrity of Output Object
test_that("Test integrity of output object", {
  data <- data.frame(aval = factor(c(1, 0)), trtp = factor(c(1, 0)), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data)
  fit1 <- predict_counterfactuals(fit, "trtp")
  original_fit <- fit1
  modified_fit <- average_predictions(fit1)
  expect_equal(length(names(modified_fit)), length(names(original_fit)) + 1) # Check for only one new addition
})

# Test Case 5: Handling Different Data Types
test_that("Test handling of different data types", {
  data <- data.frame(aval = factor(c(1, 0)), trtp = factor(c("A", "B")), bl_cov = rnorm(100))
  fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data)
  fit1 <- predict_counterfactuals(fit, "trtp")
  result <- average_predictions(fit1)
  expect_type(result$counterfactual.means, "double")
  expect_length(levels(.get_data(fit)[["trtp"]]), length(result$counterfactual.means))
})

test_that("Check model is sanitized", {
  fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
    predict_counterfactuals(trt = "trtp")
  fit1$sanitized <- NULL
  expect_warning(
    average_predictions(object = fit1),
    "Input object did not meet the expected format for beeca. Results may not be valid. Consider using ?get_marginal_effect",
    fixed = TRUE
  )
})
