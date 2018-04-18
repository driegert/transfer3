# Utility functions

#' Standardizes each column of a data.frame
#' Standardizes by subtracting the mean and dividing by the standard deviation.
#'
#' @export
standardizeDf <- function(x){
  stopifnot(class(x) == "data.frame")

  data.frame(lapply(x, function(y){ (y - mean(y)) / sd(y) }))
}

#' @export
blockStartIdx <- function(n, blockSize, overlap){
  increment <- ceiling(blockSize * (1-overlap))

  sectIdx <- seq(1, n-blockSize+1, by=increment)
  numSect <- length(sectIdx)

  list(startIdx = sectIdx, increment = increment, nBlock = numSect)
}
