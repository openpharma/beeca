# set up example ----------------------------------------------------------
data <- trial01
data$trtp <- factor(data$trtp)
data_complete <- na.omit(data)

## Fit a logistic regression working model and pass it to beeca
fit <- glm(aval ~ trtp + bl_cov, family = binomial(link = "logit"), data = data_complete) |>
  get_marginal_effect(trt="trtp", method="Ye", contrast="diff", reference = "0")


test_that("tidy_beeca stops if input is not a 'beeca' object", {
  non_beeca_object <- lm(mpg ~ wt, data = mtcars)
  expect_error(tidy_beeca(non_beeca_object), "x must be a 'beeca' object")
})


test_that("tidy_beeca returns a tibble", {
  result <- tidy_beeca(fit)
  expect_true(tibble::is_tibble(result))
})

test_that("tidy_beeca returns correct columns without confidence intervals", {
  result <- tidy_beeca(fit, conf.int = FALSE)
  expect_named(result, c("term", "contrast", "estimate", "std.error", "statistic", "p.value"))
})

test_that("tidy_beeca returns correct columns with confidence intervals", {
  result <- tidy_beeca(fit, conf.int = TRUE)
  expect_named(result, c("term", "contrast", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
})

test_that("tidy_beeca computes confidence intervals correctly", {
  result <- tidy_beeca(fit, conf.int = TRUE, conf.level = 0.95)
  conf.low <- result$estimate - qnorm(1 - 0.95 / 2) * result$std.error
  conf.high <- result$estimate + qnorm(1 - 0.95 / 2) * result$std.error
  expect_equal(result$conf.low, conf.low, tolerance = 1e-8)
  expect_equal(result$conf.high, conf.high, tolerance = 1e-8)
})


test_that("tidy_beeca returns correct values", {
  result <- tidy_beeca(fit, conf.int = FALSE)
  expect_equal(result$term, "trtp")
  expect_equal(result$estimate, fit$marginal_results$STATVAL[fit$marginal_results$STAT == "diff"])
  expect_equal(result$std.error, fit$marginal_results$STATVAL[fit$marginal_results$STAT == "diff_se"])
  expect_equal(result$statistic, NA)
  expect_equal(result$p.value, NA)
})

test_that("tidy_beeca returns correct values with confidence intervals", {
  result <- tidy_beeca(fit, conf.int = TRUE)
  expect_equal(result$term, "trtp")
  expect_equal(result$estimate, fit$marginal_results$STATVAL[fit$marginal_results$STAT == "diff"])
  expect_equal(result$std.error, fit$marginal_results$STATVAL[fit$marginal_results$STAT == "diff_se"])
  expect_equal(result$statistic, NA)
  expect_equal(result$p.value, NA)
  conf.low <- result$estimate - qnorm(1 - 0.95 / 2) * result$std.error
  conf.high <- result$estimate + qnorm(1 - 0.95 / 2) * result$std.error
  expect_equal(result$conf.low, conf.low, tolerance = 1e-8)
  expect_equal(result$conf.high, conf.high, tolerance = 1e-8)
})



