% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merging_functions.R
\name{merge}
\alias{merge}
\title{Compute p-values after merging each possible pair of X and Y labels.}
\usage{
merge(
  T_vector,
  Sigma,
  dx,
  dy,
  metric,
  xsearch,
  ysearch,
  colsample_bylevel,
  trees = list(NULL, NULL),
  categories = list(NULL, NULL)
)
}
\arguments{
\item{T_vector}{vector of length dx*dy.}

\item{Sigma}{matrix of dimension dx\emph{dy by dx}dy.}

\item{dx}{integer dimension of X.}

\item{dy}{integer dimension of Y.}

\item{metric}{str controlling the metric used for searching.}

\item{xsearch}{str "greedy", "ordinal", or "tree" controlling the X search procedure.}

\item{ysearch}{str "greedy", "ordinal", or "tree" controlling the Y search procedure.}

\item{colsample_bylevel}{float in (0,1]. The subsample ratio of potential merges to be considered.
Higher values result in more regularisation.}

\item{trees}{list containing tree structures (only used if xsearch or ysearch == tree)}

\item{categories}{list containing current merging structures (only used if xsearch or ysearch == tree)}
}
\value{
matrix with 4 named columns. "dim" equals 1 for X-label merges and 2 for Y-label merges.
ind1 < ind2 correspond to the indices of the labels being merged. "value"
is a measure of belief in the null hypothesis of conditional independence, with smaller
values favouring rejection.
}
\description{
X-merges are only considered if dx>2, similarly for Y.
}
