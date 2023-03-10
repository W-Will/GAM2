% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rscript.R
\name{SmoothingSpline}
\alias{SmoothingSpline}
\title{A implementation of Smoothing Spline function and get the coefficients}
\usage{
SmoothingSpline(X, y, lambda = NA, max_iter = 1000, tol = 1e-06)
}
\arguments{
\item{X:}{Attribute matrix (n*p)}

\item{y:}{Vector of response}

\item{lambda:}{Penalty vector on smoothness of each attribute, default to be all 0}

\item{max_iter:}{Maximum iterations}

\item{tol:}{Convergence tolerence}
}
\value{
: List of coefficients containing p coefficient vectors
}
\description{
A implementation of Smoothing Spline function and get the coefficients
}
