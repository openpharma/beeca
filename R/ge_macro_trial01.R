#' Output from the Ge et al (2011) SAS macro applied to the trial01 dataset
#'
#' For purposes of implementation comparisons, these are the result outputs from
#' the SAS macro provided with the Ge et al (2011) publication
#' (https://doi.org/10.1177/009286151104500409), applied to the trial01 dataset
#' included with beeca, adjusting for treatment (trtp) and a single covariate
#' (bl_cov) and targeting a risk difference contrast.
#'
#' @format `ge_macro_trial01`
#' A tibble with 1 row and 6 columns:
#' \describe{
#'    \item{diff}{Marginal risk difference estimate}
#'    \item{se}{Standard error of marginal risk difference estimate}
#'    \item{pt}{Marginal risk in treated}
#'    \item{pC}{Marginal risk in controls}
#'    \item{lower}{Lower bound of 95 percent confidence interval of risk difference estimate}
#'    \item{upper}{Upper bound of 95 percent confidence interval of risk difference estimate}
#' }
#'
"ge_macro_trial01"
