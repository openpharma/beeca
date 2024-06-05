# Preparing data and fitting a GLM model
trial01$trtp <- factor(trial01$trtp)
fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01[-1, ])

# Generating counterfactual predictions
fit1 <- predict_counterfactuals(fit1, "trtp")

# Predictions 1:5 from fit1$counterfactual.predictions for trial01 data
sample_df <- data.frame(
  cf_pred_0 = c(0.533, 0.537, 0.481, 0.510, 0.428),
  cf_pred_1 = c(0.463, 0.468, 0.413, 0.441, 0.362)
)

test_that("Correct counterfactual prediction for trial01 data", {
  expect_equal(round(fit1$counterfactual.predictions[1:5, ], 3)[[1]], sample_df[[1]])
  expect_equal(round(fit1$counterfactual.predictions[1:5, ], 3)[[2]], sample_df[[2]])
})
