# Binary Endpoint Estimation with Covariate Adjustment

<!-- badges: start -->

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) 
![CRAN status](https://www.r-pkg.org/badges/version/beeca)
[![R-CMD-check](https://github.com/openpharma/beeca/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/openpharma/beeca/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of **beeca** is to provide an implementation solution with a simple user interface to estimate marginal estimands in a binary endpoint setting with covariate adjustment. The primary aim of this lightweight implementation is to facilitate quick industry adoption and use within GxP environments. A secondary aim is to support the simulation studies included in the manuscript [Magirr et al. (2024)](https://osf.io/9mp58/). 


## Installation

Type | Source | Command
---|---|---
Release | CRAN | `install.packages("beeca")`
Development | GitHub | `remotes::install_github("openpharma/beeca")`

## Methodology

Motivated by the recent [FDA guidance (2023)](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/adjusting-covariates-randomized-clinical-trials-drugs-and-biological-products) on "Adjusting for Covariates in Randomized Clinical Trials for Drugs and Biological Products Guidance for Industry"" and its recommendations on robust variance estimation, we implemented two approaches, namely [Ge et al. (2011)](https://link.springer.com/article/10.1177/009286151104500409) and [Ye et al. (2023)](https://doi.org/10.1080/24754269.2023.2205802), to perform variance estimation for a treatment effect estimator based on g-computation under a logistic regression working model.

## Scope

The package is designed to estimate marginal (unconditional) estimands in a binary endpoint setting with covariate adjustment. It is suited for 2-arm clinical trials with or without covariate adaptive (stratified permuted block or biased coin) randomization where the summary measure of the marginal estimand is one of (risk difference, odds ratio, risk ratio, log odds ratio, log risk ratio). For practical considerations on the implications of covariate adjustment in superiority vs non-inferiority trials, please see [Nicholas et al. (2015)](https://doi.org/10.1002/sim.6447) and [Morris et al. (2022)](https://doi.org/10.1186/s13063-022-06097-z).

## Example

This is a basic example which shows how to obtain the point and variance estimates of a marginal estimand with covariate adjusted working model:

``` r
library(beeca)

## Set treatment to a factor
trial01$trtp <- factor(trial01$trtp)

## Fit a logistic regression working model and pass it to beeca
fit1 <- glm(aval ~ trtp + bl_cov, family="binomial", data=trial01) |>
  get_marginal_effect(trt="trtp", method="Ye", contrast="diff")

## View the results in Analysis Results Data (ARD) structure
fit1$marginal_results
```

## Package documentation 

The package documentation can be found [here](https://openpharma.github.io/beeca/). For a brief overview of the different estimands and their estimation, please see vignette [`vignette("estimand_and_implementations")`](https://openpharma.github.io/beeca/articles/estimand_and_implementations.html).

## Quality checks

Where possible we have cross checked the {beeca} package with alternative implementations in `SAS` and `R`. For example, the Ge et al. method which applies the delta method has been cross checked against the `SAS` [%margins](https://support.sas.com/kb/63/038.html) macro and the R packages [{margins}](https://cran.r-project.org/package=margins) and [{marginaleffects}](https://cran.r-project.org/package=marginaleffects). The Ye et al. method has been cross checked against [{RobinCar}](https://cran.r-project.org/package=RobinCar/). 

## Package authors

-   [Alex Przybylski](mailto:alexander.przybylski@novartis.com?subject=beeca){.email}
-   [Craig Wang](mailto:craig.wang@novartis.com?subject=beeca){.email}
-   [Dominic Magirr](mailto:dominic.magirr@novartis.com?subject=beeca){.email}
-   [Mark Baillie](mailto:mark.baillie@novartis.com?subject=beeca){.email}

## Acknowledgments

Our lightweight implementation was inspired and aided by the more comprehensive [{RobinCar}](https://cran.r-project.org/package=RobinCar/) package, developed by
Marlena Bannick, Ting Ye et al. We thank the [ASA-BIOP Covariate Adjustment Scientific Working Group](https://carswg.github.io/) for valuable feedback and discussions. 

Further development of covariate adjustment software is by the [Software Subteam](https://carswg.github.io/subteam_software.html) of ASA-BIOP Covariate Adjustment Scientific Working Group.

## References

* FDA. 2023. "Adjusting for Covariates in Randomized Clinical Trials for Drugs and Biological Products. Final Guidance for Industry."  <https://www.fda.gov/regulatory-information/search-fda-guidance-documents/adjusting-covariates-randomized-clinical-trials-drugs-and-biological-products>

* Ge, Miaomiao, L Kathryn Durham, R Daniel Meyer, Wangang Xie, and Neal Thomas. 2011. "Covariate-Adjusted Difference in Proportions from Clinical Trials Using Logistic Regression and Weighted Risk Differences." *Drug Information Journal: DIJ/Drug Information Association* 45: 481--93. <https://doi.org/10.1177/009286151104500409>

* Magirr, Dominic, Mark Baillie, Craig Wang, and Alexander Przybylski. 2024. “Estimating the Variance of Covariate-Adjusted Estimators of Average Treatment Effects in Clinical Trials with Binary Endpoints.” OSF. May 16. <https://osf.io/9mp58>.

* Ye, Ting, Marlena Bannick, Yanyao Yi, and Jun Shao. 2023. "Robust Variance Estimation for Covariate-Adjusted Unconditional Treatment Effect in Randomized Clinical Trials with Binary Outcomes." *Statistical Theory and Related Fields* 7 (2): 159--63. <https://doi.org/10.1080/24754269.2023.2205802>
