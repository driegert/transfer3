% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{blockStartIdx}
\alias{blockStartIdx}
\title{Determine block starting indices}
\usage{
blockStartIdx(n, blockSize, overlap)
}
\value{
A \code{list} containing the block starting indices, the increment between blocks,
and the number of total blocks given a vector size of n.
}
\description{
Determines starts using blockSize and overlap between blocks
}
