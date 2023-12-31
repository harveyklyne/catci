% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other_methods.R
\name{chi_sq}
\alias{chi_sq}
\title{Evaluate the chi-squared p-value for T_vector given Sigma.}
\usage{
chi_sq(T_vector, Sigma, dx, dy)
}
\arguments{
\item{T_vector}{vector of length dx*dy.}

\item{Sigma}{matrix of dimension dx\emph{dy by dx}dy.}

\item{dx}{integer dimension of X.}

\item{dy}{integer dimension of Y.}
}
\value{
float between 0 and 1
}
\description{
Under the null T is asymptotically N(0, Sigma), and so T^T \%\emph{\% Sigma^+ \%}\% T
follows a chi-square distribution with (dx-1)(dy-1) degrees of freedom, where
Sigma^+ is the pseudoinverse of Sigma.
}
