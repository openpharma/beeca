% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sanitize.R
\name{sanitize_model}
\alias{sanitize_model}
\title{(internal) Sanitize functions to check model and data within GLM model object}
\usage{
sanitize_model(model, ...)
}
\arguments{
\item{model}{a model object, currently only
\link[stats]{glm} with binomial family canonical link is supported.}

\item{...}{arguments passed to or from other methods.}
}
\value{
if model is non-compliant will throw warnings or errors.
}
\description{
Performs checks on a GLM model object to ensure it meets specific criteria
required for further analysis using other functions from the \code{beeca} package.

This includes verifying the model's family, link function, data completeness
and mode convergence.

Currently it supports models with a binomial family and canonical logit link.
}
\keyword{internal}
