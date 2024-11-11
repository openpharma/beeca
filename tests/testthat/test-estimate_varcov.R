# set up examples ----------------------------------------------------------

data01 <- trial01 |>
  transform(trtp = as.factor(trtp)) |>
  dplyr::filter(!is.na(aval))
fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
  beeca:::predict_counterfactuals(trt = "trtp") |>
  beeca:::average_predictions()

gr <- c(-1, 1)

# 3 arm example
data02 <- trial02_cdisc |>
  transform(TRTPN = as.factor(TRTPN))

fit_3arm <- glm(AVAL ~ TRTPN + SEX, family = "binomial", data = data02) |>
  beeca::predict_counterfactuals(trt = "TRTPN") |>
  beeca::average_predictions()


# test warnings/errors ------------------------------------------------------------


test_that("Correctly throwing errors on missing components", {
  expect_error(
    glm(aval ~ trtp + bl_cov, family = "binomial", data = data01) |>
      estimate_varcov(method = "Ye")
  )
})

test_that("Correctly throwing errors on argument mismatch", {
  expect_error(
    fit1 |>
      estimate_varcov(method = "Fe")
  )
})

test_that("Correctly throwing errors on missing stratification variable", {
  expect_error(
    fit1 |>
      estimate_varcov(method = "Ye", strata = "bl_cov_c")
  )
})

test_that("Correctly throwing warnings when strata supplied with method Ge", {
  expect_warning(
    fit1 |>
      estimate_varcov(method = "Ge", strata = "bl_cov_c")
  )
})

test_that("Check model is sanitized", {
  fit1$sanitized <- NULL

  expect_warning(
    estimate_varcov(object = fit1),
    "Input object did not meet the expected format for beeca. Results may not be valid. Consider using ?get_marginal_effect",
    fixed = TRUE
  )
})

# test results ------------------------------------------------------------

## Ye ----------------------------------------------------------------------

V1 <- estimate_varcov(fit1, method = "Ye")$robust_varcov

V1_3arm <- estimate_varcov(fit_3arm, method = "Ye")$robust_varcov

# if RobinCar is available
robincar_available <- requireNamespace("RobinCar", versionCheck = list(name = "RobinCar", op = "==", version = "0.3.0"), quietly = T)

if (robincar_available){
  test_that("Correct variance calculation for Ye's method matching RobinCar", {
    expect_equal(
      (t(gr) %*% V1 %*% gr)[1, 1],
      RobinCar::robincar_glm(data.frame(fit1$data),
                             response_col = "aval",
                             treat_col = "trtp",
                             formula = fit1$formula,
                             g_family = fit1$family,
                             contrast_h = "diff"
      )$contrast$varcov[1,1]
    )
  })


  test_that("3-arm: Correct variance calculation for Ye's method matching RobinCar", {
    expect_equal(
      matrix(V1_3arm, nrow=3, ncol=3),
      unname(RobinCar::robincar_glm(data.frame(fit_3arm$data),
                                    response_col = "AVAL",
                                    treat_col = "TRTPN",
                                    formula = fit_3arm$formula,
                                    g_family = fit_3arm$family)$varcov)
    )
  })


} else {
  test_that("Correct variance calculation for Ye's method matching RobinCar", {
    expect_equal(
      round((t(gr) %*% V1 %*% gr)[1, 1], 9),
      0.003706813
    )
  })

  test_that("3-arm: Correct variance calculation for Ye's method matching RobinCar", {
    expect_equal(
      matrix(V1_3arm, nrow=3, ncol=3),
      matrix(c(2.652272364e-03, 2.576978801e-06, -5.244417729e-07,
               2.576978801e-06, 2.303420439e-03,  1.053457633e-05,
               -5.244417729e-07, 1.053457633e-05,  2.379076091e-03),
             nrow=3, ncol=3)
    )
  })

}


# varcov must be symmetric
test_that("Ye varcov is symmetric", {
  expect_equal(
    V1[lower.tri(V1)],
    V1[upper.tri(V1)]
  )
})



## Ye - stratification -----------------------------------------------------

fit2 <- glm(aval ~ trtp + bl_cov + bl_cov_c2, data = data01 |> transform(bl_cov_c2 = cut(bl_cov, 4)), family = "binomial") |>
  predict_counterfactuals(trt = "trtp") |>
  average_predictions()

expect_warning(
  V2 <- estimate_varcov(fit2, method = "Ye", strata = "bl_cov_c2")$robust_varcov,
  "More than three unique values are found in stratification variable bl_cov_c2. Please double check if the correct stratification variables are supplied via `strata` argument.",
  fixed = TRUE
)

V2_3arm <- estimate_varcov(fit_3arm, method = "Ye", strata = "SEX")$robust_varcov

if (robincar_available){
  test_that("Correct variance calculation for Ye's method matcing RobinCar", {
    expect_equal(
      (t(gr) %*% V2 %*% gr)[1, 1],
      RobinCar::robincar_glm(data.frame(fit2$data),
                             response_col = as.character(fit2$formula[2]),
                             treat_col = "trtp",
                             formula = fit2$formula,
                             car_strata_cols = "bl_cov_c2",
                             car_scheme = "permuted-block",
                             g_family = fit2$family,
                             contrast_h = "diff",
      )$contrast$varcov[1,1]
    )
  })

  test_that("3-arm: Correct variance calculation for Ye-stratified method matching RobinCar", {
    expect_equal(
      matrix(V2_3arm, nrow=3, ncol=3),
      unname(RobinCar::robincar_glm(data.frame(fit_3arm$data),
                                    response_col = "AVAL",
                                    treat_col = "TRTPN",
                                    formula = fit_3arm$formula,
                                    car_strata_cols = "SEX",
                                    car_scheme = "permuted-block",
                                    g_family = fit_3arm$family)$varcov)
    )
  })

} else {
  test_that("Correct variance calculation for Ye's method matching RobinCar", {
    expect_equal(
      round((t(gr) %*% V2 %*% gr)[1, 1], 9),
      0.003664811
    )
  })


  test_that("3-arm: Correct variance calculation for Ye-stratified method matching RobinCar", {
    expect_equal(
      matrix(V2_3arm, nrow=3, ncol=3),
      matrix(c(2.622886513e-03, -7.797320489e-06, -4.923527008e-06,
               -7.797320489e-06,  2.288915038e-03,  1.360521021e-05,
               -4.923527008e-06,  1.360521021e-05,  2.376256827e-03),
             nrow=3, ncol=3)
    )
  })

}

# varcov must be symmetric
test_that("Ye-stratified varcov is symmetric", {
  expect_equal(
    V2[lower.tri(V2)],
    V2[upper.tri(V2)]
  )
})



## Ye - mod ----------------------------------------------------------------

# need a smaller dataset to obtain different results with mod argument
fit3 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = data01[1:150, ]) |>
  beeca:::predict_counterfactuals(trt = "trtp") |>
  beeca:::average_predictions()

V3 <- estimate_varcov(fit3, method = "Ye", mod = T)$robust_varcov

test_that("Correct variance calculation for Ye's method based on original manuscript reference matching modified RobinCar", {
  expect_equal(
    round((t(gr) %*% V3 %*% gr)[1, 1], 9),
    0.006669802
    # based on function from Variance_Ye vignette: ye_var_robincar_mod(fit3, "trtp", convert_back = F)
  )
})



# varcov must be symmetric
test_that("Ye-mod varcov is symmetric", {
  expect_equal(
    V3[lower.tri(V3)],
    V3[upper.tri(V3)]
  )
})


## Ge ----------------------------------------------------------------------

ge_var_paper <- function(glmfit, trt) {
  if (is.factor(glmfit$model[, trt]) & nlevels(glmfit$model[, trt]) != 2) {
    stop("Only 2 levels of trt supported.")
  }

  pder <- function(ahat, vc, x) {
    ### ahat: logistic regression parameters
    ### vc: variance-covariance matrix of ahat
    ### x: full model matrix of the logistic regression
    ### return mean of phat on x and its se
    phat <- plogis(x %*% ahat)
    pbar <- mean(phat)
    pderiv <- t(phat * (1 - phat)) %*% x / nrow(x)
    sepbar <- sqrt(pderiv %*% vc %*% t(pderiv))
    return(list(pbar = pbar, sepbar = sepbar, pderiv = pderiv))
  }

  difP <- function(glmfit) {
    ### estimate the proportion difference and its standard error
    df <- glmfit$model

    vc <- vcov(glmfit)
    df[, trt] <- 1
    mat <- model.matrix(glmfit$formula, data = df)
    pderT <- pder(coef(glmfit), vc, mat)
    df[, trt] <- 0
    mat <- model.matrix(glmfit$formula, data = df)
    pderC <- pder(coef(glmfit), vc, mat)

    difb <- pderT$pbar - pderC$pbar
    sedif <- sqrt((pderT$pderiv - pderC$pderiv) %*% vc %*% t(pderT$pderiv - pderC$pderiv))

    return(list(
      pT = pderT$pbar,
      pC = pderC$pbar,
      dif = difb,
      sedif = sedif,
      var = sedif**2
    ))
  }

  return(difP(glmfit)$var[[1]])
}

V4 <- estimate_varcov(fit1, method = "Ge", type = "model-based")$robust_varcov

test_that("Correct variance calculation for Ge's method matching manuscript code", {
  expect_equal(
    (t(gr) %*% V4 %*% gr)[1, 1],
    ge_var_paper(fit1, "trtp")
  )
})


# varcov must be symmetric
test_that("Ge varcov is symmetric", {
  expect_equal(
    V4[lower.tri(V4)],
    V4[upper.tri(V4)]
  )
})
