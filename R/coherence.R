# coherence based methods:

#' @details Only calculates the forward coherence currently.
#' @export
coherence <- function(x, y, method = c("aveOfMsc", "mscOfAve", "mscOfAveCoh")
                      , nw = 4, k = 7, nFFTx = NULL, centre = "none"
                      , dtx = 1, dty = 1, blockSizex = NULL, overlap = 0
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
    nFFTx <- 2^ceiling(log2(nx + 1))
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
                             , convertMscToNormal = convertMscToNormal)
         , cFreq = freqx[freqRangeIdx[1]:freqRangeIdx[2]]
         , oFreq = oFreq
         , info = info)
  }
}

#' @export
aveOfMscBlock <- function(x, y, blockx, blocky
                          , freqRangeIdx, nFreqRangeIdx
                          , maxFreqOffsetIdx, nOffsets
                          , nw, k, nFFTx, nFFTy, centre, dtx, dty
                          , convertMscToNormal){
  msc <- matrix(0, nrow = nOffsets, ncol = nFreqRangeIdx)

  for (i in 1:blockx$nBlock){
    ykx <- eigenCoef(x[ blockx$startIdx[i]:(i*blockx$blockSize) ], nw = nw, k = k, nFFT = nFFTx
                     , centre = centre, deltat = dtx)
    yky <- eigenCoef(y[ blocky$startIdx[i]:(i*blocky$blockSize) ], nw = nw, k = k, nFFT = nFFTy
                     , centre = centre, deltat = dty)

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

  msc
}

# the assumption here is that the df's are the same in yk1 and yk2, which means nFFT1 = a*nFFT2 (a = dt1/dt2 an integer)
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
    out <- .Fortran("cohMsc"
                    , coh = double(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2], ])
                    , yk2 = as.complex(yk2[ (freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]+maxFreqOffsetIdx), ])
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(nrow(yk2.f))
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  }

  matrix(out$coh, nrow = nOffsets, ncol = nfreq1)
}

mscOfAve <- function(){
  -1
}

mscOfAveCoh <- function(){

}
