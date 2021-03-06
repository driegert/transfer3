% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilsSpec.R
\name{eigenCoefWithNegYk}
\alias{eigenCoefWithNegYk}
\title{Extends positive frequency eigencoefficients to include negative frequencies}
\usage{
eigenCoefWithNegYk(yk2, freqRangeIdx, maxFreqOffsetIdx, multPred = FALSE)
}
\description{
Adds conjugate eigencoefficient values to the "top" (front?) of the eigencoefficient matrix.
}
\details{
The multPred argument is assuming that you are passing in a 4D array:
freqs x tapers x blocks x predictors
}
