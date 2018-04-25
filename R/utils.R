# Utility functions

#' Standardizes each column of a data.frame
#' Standardizes by subtracting the mean and dividing by the standard deviation.
#'
#' @export
standardizeDf <- function(x){
  stopifnot(class(x) == "data.frame")

  data.frame(lapply(x, function(y){ (y - mean(y)) / sd(y) }))
}


#' Determine block starting indices
#'
#' Determines starts using blockSize and overlap between blocks
#'
#' @return A \code{list} containing the block starting indices, the increment between blocks,
#' and the number of total blocks given a vector size of n.
#'
#' @export
blockStartIdx <- function(n, blockSize, overlap){
  increment <- ceiling(blockSize * (1-overlap))

  sectIdx <- seq(1, n-blockSize+1, by=increment)
  numSect <- length(sectIdx)

  list(startIdx = sectIdx, increment = increment, nBlock = numSect, blockSize = blockSize)
}

#' Convert the MSC to Standard Normal
#'
#' @param msc A \code{matrix} containing MSC values.
#' @param k An \code{integer} representing the number of tapers used in MSC calculation.
#' @param J An \code{integer} indicating the number of blocks used in MSC calculation.
#' @param mu The mean of the Standard Normal.
#' @param sd The standard deviation of the Standard Normal.
#' @param avoidInf A \code{logical} indicating whether to avoid potential infinite values.
#' Returns 10 or -10 instead.
#'
#' Do NOT change this parameter.  Results are not reliable if changed.
#'
#' @export
msc2norm <- function(msc, k, J = 1, mu = 0, sd = 1, avoidInf = FALSE){

  tran <- 1 - (1-msc)^(J*(k-1))

  if (avoidInf && tran == 1){
    return(10)
  } else if (avoidInf && tran == -1){
    return(-10)
  } else if (tran == 1 || tran == -1){
    warning("msc could possibly cause +/- Inf value.")
  }

  mu + sqrt(2)*sd*erfinv(2*tran - 1)
}

#' @references Note used: Abromowitz & Stegun - 7.1.26 - http://people.math.sfu.ca/~cbm/aands/page_299.htm
erfinv <- function(x){

  qnorm((x + 1)/2)/sqrt(2) # this is from qnorm R documention - bottom examples.

# These do NOT seem correct - at least, they don't match erf.inv() in msc2norm() ... ^^
  # p <- 0.3275911
  # a <- c(0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429)
  # t <- 1 / (1 + p*x)
  #
  # 1 - (a[1]*t + a[2]*t^2 + a[3]*t^3 + a[4]*t^4 + a[5]*t^5)*exp(-x^2)
}
