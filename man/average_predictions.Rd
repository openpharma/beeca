% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/average_predictions.R
\name{average_predictions}
\alias{average_predictions}
\title{Average over counterfactual predictions}
\usage{
average_predictions(object)
}
\arguments{
\item{object}{a fitted \code{\link[stats]{glm}} object augmented with
counterfactual predictions named: \code{counterfactual.predictions}}
}
\value{
an updated \code{glm} object appended with an additional component
\code{counterfactual.means}.
}
\description{
\code{average_predictions()} averages counterfactual predictions stored within
a \code{glm} object. This is pivotal for estimating treatment contrasts and
associated variance estimates using g-computation. The function assumes
predictions are generated via \code{predict_counterfactuals()}.
}
\details{
The \code{average_predictions()} function calculates the average over the
counterfactual predictions which can then be used to estimate a treatment
contrast and associated variance estimate.

The function appends a \code{glm} object with the
averaged counterfactual predictions.

Note: Ensure that the \code{glm} object has been adequately prepared with
\code{predict_counterfactuals()} before applying \code{average_predictions()}.
Failure to do so may result in errors indicating missing components.
}
\examples{

# Use the trial01 dataset
data(trial01)

# ensure the treatment indicator is a factor
trial01$trtp <- factor(trial01$trtp)

# fit glm model for trial data
fit1 <- glm(aval ~ trtp + bl_cov, family = "binomial", data = trial01)

# Preprocess fit1 as required by average_predictions
fit2 <- fit1 |>
  predict_counterfactuals(trt = "trtp")

# average over the counterfactual predictions
fit3 <- average_predictions(fit2)

# display the average predictions
fit3$counterfactual.means

}
\seealso{
\code{\link[=predict_counterfactuals]{predict_counterfactuals()}} for generating counterfactual
predictions.

\code{\link[=estimate_varcov]{estimate_varcov()}} for estimating the variance-covariate matrix
of mariginal effects

\code{\link[=get_marginal_effect]{get_marginal_effect()}} for estimating marginal effects directly
from an original \code{\link[stats]{glm}} object
}
