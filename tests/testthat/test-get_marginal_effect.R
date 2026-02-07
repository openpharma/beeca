# set up example ----------------------------------------------------------

data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))

fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
  get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

# 3-arm example
data02 <- trial02_cdisc |>
  transform(TRTPN = as.factor(TRTPN))

fit_3arm <- glm(AVAL ~ TRTPN + SEX, family = "binomial", data = data02) |>
  get_marginal_effect(trt = "TRTPN", method = "Ye", contrast = "diff", reference = "1")


# test warnings/errors ------------------------------------------------------------

test_that("Correctly throwing errors on missing argument", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      get_marginal_effect(method = "Ye"))
})

test_that("Correctly throwing errors on argument mismatch", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      get_marginal_effect(trt = "trtp", method = "Fe"))
})

test_that("Correctly throwing warnings on missing value", {
  expect_warning(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01 |>
          transform(trtp = as.factor(trtp))) |>
      get_marginal_effect(trt = "trtp", method = "Ye", reference = "0"),
    "There is 1 record omitted from the original data due to missing values, please check if they should be imputed prior to model fitting."
  )
})


# test results ------------------------------------------------------------

test_that("Correctly producing results in ARD format", {
  expect_equal(dim(fit1$marginal_results), c(12, 8))
})

test_that("Correctly producing results", {
  expect_false(any(is.na(fit1$marginal_results)))
})


test_that("3-arm: Correctly producing results in ARD format", {
  expect_equal(dim(fit_3arm$marginal_results), c(19, 8))
})

test_that("3-arm: Correctly producing results", {
  expect_false(any(is.na(fit_3arm$marginal_results)))
})
