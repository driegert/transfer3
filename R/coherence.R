# coherence based methods:

#' Calculates the coherence between two series.
#'
#' Multitaper coherence estimate.
#'
#' @param x slower sampled series - the response (stupid)
#' @param y faster sampled series - the predictor
#' @param method only "aveOfMsc" and "offTaperMsc" are currently implemented
#' @param nw blap
#' @param k blap
#' @param centre The time series is centred using one of three methods:
#' expansion onto discrete prolate spheroidal sequences ('Slepian'),
#' arithmetic mean ('arithMean'), trimmed mean ('trimMean'), or not at all ('none').
#'
#'
#' @details Only calculates the forward coherence currently.  'x' is the response if
#' the series do not have the same sampling rate ... which is confusing.
#' @export
coherence <- function(x, y, method = c("aveOfMsc", "offTaperMsc", "mscOfAve", "mscOfAveCoh")
                      , nw = 4, k = 7, nFFTx = NULL, centre = "none"
                      , dtx = 1, dty = 1, adaptiveWeighting = TRUE
                      , blockSizex = NULL, overlap = 0
                      , freqRange = NULL, maxFreqOffset = 0
                      , convertMscToNormal = FALSE
                      , namex = "x", namey = "y"){
  # add some checks for NA's:
  dtRatio <- dtx / dty
  nx <- length(x)
  ny <- length(y)

  if (is.null(blockSizex)){
    blockSizex <- nx
  }
  blockSizey <- dtRatio * blockSizex

  # get correct blocksizes:
  stopifnot(dtRatio >= 1 # series y should be sampled as or more quickly than series x
            , blockSizex <= length(x) # blockSize can't be larger than the length of the series
            , dtRatio - trunc(dtRatio) == 0) # need the ratio to be an integer for this to work properly

  if (is.null(nFFTx)){
    nFFTx <- 2^ceiling(log2(nx) + 1)
  }
  nFFTy <- dtRatio * nFFTx

  freqx <- posFreq(nFFTx, dtx)
  freqy <- posFreq(nFFTy, dty)

  df <- 1/(nFFTx * dtx)

  nfreqx <- length(freqx)
  nfreqy <- length(freqy)

  if (is.null(freqRange)){
    freqRange <- c(1, 1/(2*dtx))
    freqRangeIdx <- c(1, nfreqx)
  } else {
    freqRangeIdx <- c(max(1, floor(freqRange[1] / df)), min(floor(freqRange[2] / df), nfreqx))
  }

  nFreqRangeIdx <- freqRangeIdx[2] - freqRangeIdx[1] + 1

  maxFreqOffsetIdx <- ceiling(maxFreqOffset / df)
  nOffsets <- 2*maxFreqOffsetIdx + 1
  if (nOffsets > 1){
    oFreq <- c(rev(-1*freqx[2:(maxFreqOffsetIdx+1)]), 0, freqx[2:(maxFreqOffsetIdx+1)])
  } else {
    oFreq <- 0
  }

  blockx <- blockStartIdx(n = nx, blockSize = blockSizex, overlap = overlap)
  blocky <- blockStartIdx(n = ny, blockSize = blockSizey, overlap = overlap)

  info <- list(namex = namex, namey = namey
               ,method = method, nw = nw, k = k, nFFTx = nFFTx, nFFTy = nFFTy
               , centre = centre, dtx = dtx, dty = dty, dtRatio = dtRatio
               , adaptiveWeighting = adaptiveWeighting
               , blockSizex = blockSizex, blockSizey = blockSizey
               , overlap = overlap
               , nBlocks = blockx$nBlock
               , freqRange = freqRange, freqRangeIdx = freqRangeIdx
               , nFreqRangeIdx = nFreqRangeIdx
               , maxFreqOffset = maxFreqOffset, maxFreqOffsetIdx = maxFreqOffsetIdx
               , convertMscToNormal = convertMscToNormal)

  if (method[1] == "aveOfMsc"){
    list(coh = aveOfMscBlock(x = x, y = y, blockx = blockx, blocky = blocky
                             , freqRangeIdx = freqRangeIdx, nFreqRangeIdx = nFreqRangeIdx
                             , maxFreqOffsetIdx = maxFreqOffsetIdx, nOffsets = nOffsets
                             , nw = nw, k = k
                             , nFFTx = nFFTx, nFFTy = nFFTy
                             , centre = centre, dtx = dtx, dty = dty
                             , adaptiveWeighting = adaptiveWeighting
                             , convertMscToNormal = convertMscToNormal)
         , cFreq = freqx[freqRangeIdx[1]:freqRangeIdx[2]]
         , oFreq = oFreq
         , info = info)
  } else if (method[1] == "offTaperMsc"){
    warning("Only full block will be used.")

    ykx <- eigenCoef(x, nw = nw, k = k, nFFT = nFFTx
                     , centre = centre, deltat = dtx
                     , adaptiveWeighting = adaptiveWeighting)$yk
    yky <- eigenCoef(y, nw = nw, k = k, nFFT = nFFTy
                     , centre = centre, deltat = dty
                     , adaptiveWeighting = adaptiveWeighting)$yk
    list(coh = offTaperHelper(ykx, yky, freqRangeIdx, maxFreqOffsetIdx)
         , cFreq = freqx[freqRangeIdx[1]:freqRangeIdx[2]]
         , oFreq = oFreq
         , info = info)
  }
}

#' @export
offTaperHelper <- function(yk1, yk2, freqRangeIdx, maxFreqOffsetIdx){
  if (freqRangeIdx[2] + maxFreqOffsetIdx > nrow(yk2)){
    stop("Not enough indices to accommodate the highest offset from the largest centre frequency.")
  }

  nfreq1 <- freqRangeIdx[2] - freqRangeIdx[1] + 1
  nOffsets <- 2*maxFreqOffsetIdx + 1

  if (freqRangeIdx[1] - maxFreqOffsetIdx <= 0){
    yk2.f <- eigenCoefWithNegYk(yk2 = yk2, freqRangeIdx = freqRangeIdx
                                , maxFreqOffsetIdx = maxFreqOffsetIdx)
    out <- .Fortran("offTaperMsc"
                    , coh = double(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2], ])
                    , yk2 = as.complex(yk2.f)
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(nrow(yk2.f))
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  } else {
    numRows <- length((freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]+maxFreqOffsetIdx))
    out <- .Fortran("offTaperMsc"
                    , coh = double(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2], ])
                    , yk2 = as.complex(yk2[ (freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]+maxFreqOffsetIdx), ])
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(numRows)
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  }

  matrix(out$coh, nrow = nOffsets, ncol = nfreq1)
}

#' @export
aveOfMscBlock <- function(x, y, blockx, blocky
                          , freqRangeIdx, nFreqRangeIdx
                          , maxFreqOffsetIdx, nOffsets
                          , nw, k, nFFTx, nFFTy, centre, dtx, dty
                          , adaptiveWeighting, convertMscToNormal){

  msc <- matrix(0, nrow = nOffsets, ncol = nFreqRangeIdx)

  for (i in 1:blockx$nBlock){
    ykx <- eigenCoef(x[ blockx$startIdx[i]:(blockx$startIdx[i] + blockx$blockSize-1) ]
                     , nw = nw, k = k, nFFT = nFFTx
                     , centre = centre, deltat = dtx
                     , adaptiveWeighting = adaptiveWeighting)$yk
    yky <- eigenCoef(y[ blocky$startIdx[i]:(blocky$startIdx[i] + blocky$blockSize-1) ]
                     , nw = nw, k = k, nFFT = nFFTy
                     , centre = centre, deltat = dty
                     , adaptiveWeighting = adaptiveWeighting)$yk

    if (convertMscToNormal){
      msc <- msc + msc2norm(aveOfMsc(yk1 = ykx, yk2 = yky
                                     , freqRangeIdx = freqRangeIdx
                                     , maxFreqOffsetIdx = maxFreqOffsetIdx)
                            , k = k, avoidInf = TRUE)
    } else {
      msc <- msc + aveOfMsc(yk1 = ykx, yk2 = yky
                            , freqRangeIdx = freqRangeIdx
                            , maxFreqOffsetIdx = maxFreqOffsetIdx)
    }
  }

  msc / blockx$nBlock
}

# the assumption here is that the df's are the same in yk1 and yk2,
# which means nFFT1 = a*nFFT2 (a = dt1/dt2 an integer)

## yk1 is the response
## yk2 is the "predictor" (faster sampling rate)
aveOfMsc <- function(yk1, yk2, freqRangeIdx, maxFreqOffsetIdx){
  if (freqRangeIdx[2] + maxFreqOffsetIdx > nrow(yk2)){
    stop("Not enough indices to accommodate the highest offset from the largest centre frequency.")
  }

  nfreq1 <- freqRangeIdx[2] - freqRangeIdx[1] + 1
  nOffsets <- 2*maxFreqOffsetIdx + 1

  if (freqRangeIdx[1] - maxFreqOffsetIdx <= 0){
    yk2.f <- eigenCoefWithNegYk(yk2 = yk2, freqRangeIdx = freqRangeIdx
                                , maxFreqOffsetIdx = maxFreqOffsetIdx)
    out <- .Fortran("cohMsc"
                    , coh = double(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2], ])
                    , yk2 = as.complex(yk2.f)
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(nrow(yk2.f))
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  } else {
    numRows <- length((freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]+maxFreqOffsetIdx))
    out <- .Fortran("cohMsc"
                    , coh = double(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2], ])
                    , yk2 = as.complex(yk2[ (freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]+maxFreqOffsetIdx), ])
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(numRows)
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  }

  matrix(out$coh, nrow = nOffsets, ncol = nfreq1)
}

#' used in the transfer function functions
#' @export
mscFromEigenHelper <- function(yk1, yk2, k, nOffsets){
  out <- .Fortran("cohMsc"
                  , coh = double(nOffsets * nrow(yk1))
                  , yk1 = as.complex(yk1)
                  , yk2 = as.complex(yk2)
                  , nfreq1 = as.integer(nrow(yk1))
                  , nfreq2 = as.integer(nrow(yk2))
                  , k = as.integer(k)
                  , nOffsets = as.integer(nOffsets))

  # convert to standard normal:
  matrix(msc2norm(out$coh, k = k, avoidInf = TRUE), nrow = nOffsets, ncol = nrow(yk1))
}

mscOfAve <- function(){
  -1
}

mscOfAveCoh <- function(){

}
