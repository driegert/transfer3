#' Calculate weighted eigencoefficients
#'
#' Calculates weighted eigencoefficients for a single vector / block
#'
#' @param x A vector containing a time-series for which to calculat the eigencoefficnets.
#'
#' @details I think there may be some overhead that could be removed in the spec.mtm code.
#' But also, perhaps not...
#'
#' @export

eigenCoef <- function(x, nw = 4, k = 7, nFFT = "default", centre = "none"
                      , deltat = 1, dtUnits = "second"
                      , adaptiveWeighting = TRUE, dpssIN = NULL){

  if (adaptiveWeighting){
    tmp <- spec.mtm(x, nw = nw, k = k, nFFT = nFFT, centre = centre
                    , dpssIN = dpssIN, deltat = deltat, dtUnits = dtUnits
                    , adaptiveWeighting = adaptiveWeighting
                    , plot = FALSE, returnInternals = TRUE)$mtm[ c("eigenCoefs", "eigenCoefWt") ]
    yk <- tmp$eigenCoefs * sqrt(tmp$eigenCoefWt)
  } else {
    yk <- spec.mtm(x, nw = nw, k = k, nFFT = nFFT, centre = centre
                          , dpssIN = dpssIN, deltat = deltat, dtUnits = dtUnits
                          , adaptiveWeighting = adaptiveWeighting, Ftest = TRUE
                          , plot = FALSE, returnInternals = TRUE)$mtm$eigenCoefs
  }

  yk
}

#' @export
eigenCoefFit <- function(H, Hinfo, ykx, prednames, centFreqIdx){
  # browser()
  fitted <- array(0, dim = c(length(centFreqIdx), dim(ykx)[2], dim(ykx)[3], 1))

  for (hrow in 1:dim(Hinfo)[1]){ # For each unique offset
    for (b in 1:dim(ykx)[3]){ # by block
      fitted[, , b, 1] <- fitted[, , b, 1] + apply(ykx[centFreqIdx + Hinfo$idxOffset[hrow]
                                                       , #all columns
                                                       , b
                                                       , which(prednames == Hinfo$predictor[hrow])
                                                       , drop = FALSE]
                                                   , MARGIN = 2
                                                   , FUN = "*"
                                                   , H[, Hinfo$Hcolumn[hrow] ])
      # + mapply("*", ykx[centFreqIdx + Hinfo$idxOffset[hrow]
      #                                                                   ,
      #                                                                   , b
      #                                                                   , which(prednames == Hinfo$predictor[hrow])
      #                                                                   , drop = FALSE]
      #                                                          , H[, Hinfo$Hcolumn[hrow] ])
    }
  }

  fitted
}


#' @export
posFreq <- function(nFFT, dt){
  seq(0, 1/(2*dt), by = 1/(dt*nFFT))
}

#' Extends positive frequency eigencoefficients to include negative frequencies
#'
#' Adds conjugate eigencoefficient values to the "top" (front?) of the eigencoefficient matrix.
#'
#' @details The multPred argument is assuming that you are passing in a 4D array:
#' freqs x tapers x blocks x predictors
#'
#' @export
eigenCoefWithNegYk <- function(yk2, freqRangeIdx, maxFreqOffsetIdx, multPred = FALSE){
  if(multPred){
    abind(Conj(yk2[rev(2:(maxFreqOffsetIdx - freqRangeIdx[1] + 2)), , , , drop = FALSE])
          , yk2[1:(freqRangeIdx[2] + maxFreqOffsetIdx), , , , drop = FALSE], along = 1)
  } else {
    rbind(Conj(yk2[rev(2:(maxFreqOffsetIdx - freqRangeIdx[1] + 2)), ])
          , yk2[1:(freqRangeIdx[2] + maxFreqOffsetIdx), ])
  }
}

#' this only works for a very specific situation
#' @export
fixDeadBand <- function(msc, zeroOffsetIdx, freqRangeIdx, band, df, replaceWith = 0){
  # if ((band[1] > 0 | band[2] < 0) | (freqRangeIdx[1] != 1)){
  #   stop("fixDeadBand doesn't work unless you start at 0 and you're taking a band around 0.")
  # }

  ## assuming symmetry around f = 0:
  nIdx <- ceiling(band / df)

  top <- zeroOffsetIdx + nIdx
  bottom <- zeroOffsetIdx - nIdx

  # everything is backwards - we need the reverse diagonals.  Also, in fields::image.plot(), the
  # y-axis is "reversed" relative to matrix indices as it were...
  # so thinking about this based on those plots is difficult..
  matTop <- lower.tri(matrix(0, top, top), diag = TRUE)
  matBot <- upper.tri(matrix(0, nrow = bottom, ncol = top), diag = FALSE)
  matTop[(top-bottom+1):top, 1:top] <- matBot & matTop[(top-bottom+1):top, 1:top]

  msc[top:1, 1:top][matTop] <- replaceWith

  msc
}


convolveFilter <- function(filter, series, sides = 2){

}

#' Plots the impulse responses to a PDF
#'
#' Mostly used when using many offsets.
#'
#' @param ir A \code{matrix} with each column containing the impusles response
#' @param filename A \code{character} string indicating where to store the pdf.  Should have a
#' .pdf extension.
#' @param hTrim An \code{integer} indicating the length of each side of the filter if sides = 2
#' (i.e., L = 2*hTrim + 1), or the number of filter coefficients if sides = 1.
#' @param sides An \code{integer} indicating whether this is a one- or two-sided filter.
#'
#' @details TODO: Maybe add the actual offset label to the plots.
#'
#' @export

plotIr <- function(ir, filename, hTrim = NULL, sides = 2){
  if (is.null(hTrim)){
    hTrim <- floor(dim(ir)[1] / 2)
  }

  xlabel <- seq(-hTrim, hTrim)

  pdf(filename, height = 6, width = 8)
  par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))

  for (i in 1:ncol(ir)){
    plot(xlabel, ir[, i], type = 'l', xlab = "Lag", ylab = "Impulse Response")
  }

  dev.off()
}

#' Apply a linear filter to a time-series
#'
#' Assumes a two-sided filter.
#'
#' @param x A \code{vector} containing the series to filter.
#' @param filter A \code{vector} containing the filter coefficients.  See details.
#'
#' @details Ends up being twice as fast as stats::filter().  Only 2-sided.  Only convolution.
#'
#' The argument, filter, should be of length 2*L + 1 (i.e., odd).  The positive lag (causal) filter
#' coefficients should be in elements (L+1):(2*L) with the negative lag (acausal) filter coefficients
#' contained in the 1:L elements of filter.  The L+1 index of filter is lag-0.  This is the same convention as
#' required by \code{stats::filter()}.
#'
#' @export
filter.twoSided <- function(x, filter){
  stopifnot(length(filter) %% 2 == 1)

  out <- .Fortran("twoSidedFilter"
                  , filter = as.double(filter)
                  , series = as.double(x)
                  , seriesOut = double(length(x))
                  , hTrim = as.integer(floor(length(x) / 2))
                  , n = as.integer(length(x)))$seriesOut

  out[which(out == -999.99)] <- NA

  out
}

# DJT Paper - SPIE 1993 - Nonstationary fluctuations in stationary data
#' @export
invertEigenCoef <- function(yk, v, dt, nFFT, N){
  nfreqs <- nFFT/2+1
  cft <- rbind(yk, Conj(yk[(nfreqs-1):2, ]))
  # invert
  inv <- apply(cft, MAR = 2
               , FUN = function(x) { 1 / nFFT * Re(fft(x, inverse = TRUE))[1:N] })

  #  Invert by using the trick of the Ftest: \sum_{k} x * v_k^2 / \sum_{k} v_k^2
  scaleFactor <- v * sqrt(dt)
  fixedX2 <- rowSums(inv * scaleFactor) / rowSums(scaleFactor * scaleFactor)

  fixedX2


  # original Wes code sent through gchat:
  # cft <- rbind(egnC, Conj(egnC[(s2$mtm$nfreqs-1):2, ]))
  # # invert
  # inv <- apply(cft, MAR = 2
  #              , FUN = function(x) { 1 / s2$mtm$nFFT * Re(fft(x, inverse = TRUE))[1:N] })
  #
  # #  Invert by using the trick of the Ftest: \sum_{k} x * v_k^2 / \sum_{k} v_k^2
  # scaleFactor <- s2$mtm$dpss$v * sqrt(s2$mtm$deltaT)
  # fixedX2 <- rowSums(inv * scaleFactor) / rowSums(scaleFactor * scaleFactor)
}
