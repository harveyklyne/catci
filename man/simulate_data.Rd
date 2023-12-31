% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_settings.R
\name{simulate_data}
\alias{simulate_data}
\title{Generate an i.i.d. (X,Y,Z) dataset.}
\usage{
simulate_data(
  n,
  xnum_class,
  ynum_class,
  xsetting,
  ysetting,
  strength,
  intsetting,
  permute
)
}
\arguments{
\item{n}{integer number of samples.}

\item{xnum_class}{integer number of classes for discrete variable X.}

\item{ynum_class}{integer number of classes for discrete variable Y.}

\item{xsetting}{string determining X | Z setting.}

\item{ysetting}{string determining Y | Z setting.}

\item{strength}{float between 0 and 1 controlling the conditional dependence of Y on X given Z.}

\item{intsetting}{string determining conditional dependence setting.}

\item{permute}{logical determining if X, Y labels should be randomly permuted (violating ordinal structure)}
}
\value{
list containing (X, Y, Z) samples and potentially also corresponding population (f(Z), g(Z)).
}
\description{
Generate an i.i.d. (X,Y,Z) dataset.
}
