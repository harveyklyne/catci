% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regression_functions.R
\name{fit_xgboost}
\alias{fit_xgboost}
\title{Fit pre-tuned XGBoost regression for use in simulations.}
\usage{
fit_xgboost(X, y, num_class = NULL, params)
}
\arguments{
\item{X}{matrix of covariates.}

\item{y}{integer vector of responses in {1, ..., num_class}.}

\item{num_class}{integer number of levels of desired output.}

\item{params}{XGBoost hyperparameters.}
}
\value{
list containing a function "fit" which takes matrix input of the
same width as X, and returns a matrix of probability estimates, with each
column corresponding to a label of y.
}
\description{
Fit pre-tuned XGBoost regression for use in simulations.
}
