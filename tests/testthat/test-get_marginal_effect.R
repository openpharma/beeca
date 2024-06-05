# set up example ----------------------------------------------------------

data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))

fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
  get_marginal_effect(trt = "trtp", method = "Ye", contrast = "diff", reference = "0")

# test warnings/errors ------------------------------------------------------------

test_that("Correctly throwing errors on missing argument", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      get_marginal_effect(method = "Ye"), 'argument "trt" is missing, with no default'
  )
})

test_that("Correctly throwing errors on arugment mismatch", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      get_marginal_effect(method = "Fe"), 'argument "trt" is missing, with no default'
  )
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
