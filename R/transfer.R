# transfer function related code:

#' Calculates some transfer functions.
#'
#' @param y A single column \code{data.frame} containing the response time-series.
#' @param x A \code{data.frame} containing the predictor time-series.
#' @param decorrelate A vector containing the names of the predictors (columns of x)
#' in the order of the iterative regression.
#' @param nOffAllowed An \code{integer} representing the number of offsets, \emph{including
#' the zero-offset}, to be used in the transfer function estimation.
#' @export
tf <- function(y, x, nw = 4, k = 7, nFFTy = NULL, centre = "none"
               , dty = 1, dtx = 1, blockSizey = NULL, overlap = 0
               , cohSigLev = 0.99, nOffAllowed = 1
               , forceZeroOffset = TRUE
               , freqRange = NULL, maxFreqOffset = 0, deadBand = NULL
               , decorrelate = NULL){
  ## check all the parameters here...
  # add some checks for NA's here
  #################
  #
  #
  #
  #################

  if (forceZeroOffset & nOffAllowed == 1 & maxFreqOffset > 0){
    warning("You are requiring the zero-offset and a maximum of 1 offset (i.e., the zero-offset). \n
            Setting maxFreqOffset to 0.")
    maxFreqOffset <- 0
  }

  dtRatio <- dty / dtx
  nx <- nrow(x)
  ny <- nrow(y)

  if (is.null(blockSizey)){
    blockSizey <- ny
  }
  blockSizex <- dtRatio * blockSizey

  # get correct blocksizes:
  stopifnot(dtRatio >= 1 # series y should be sampled as or more quickly than series x
            , blockSizey <= ny # blockSize can't be larger than the length of the series
            , dtRatio - trunc(dtRatio) == 0) # need the ratio to be an integer for this to work properly

  if (is.null(nFFTy)){
    nFFTy <- 2^ceiling(log2(ny + 1))
  }
  nFFTx <- dtRatio * nFFTy

  freqy <- posFreq(nFFTy, dty)
  freqx <- posFreq(nFFTx, dtx)

  df <- 1/(nFFTy * dty)

  nfreqy <- length(freqy)
  nfreqx <- length(freqx)

  if (is.null(freqRange)){
    freqRange <- c(1, 1/(2*dty))
    freqRangeIdx <- c(1, nfreqy)
  } else {
    freqRangeIdx <- c(max(1, floor(freqRange[1] / df)), min(floor(freqRange[2] / df), nfreqy))
  }

  nFreqRangeIdx <- freqRangeIdx[2] - freqRangeIdx[1] + 1

  maxFreqOffsetIdx <- ceiling(maxFreqOffset / df)
  nOffsets <- 2*maxFreqOffsetIdx + 1
  if (nOffsets > 1){
    oFreq <- c(rev(-1*freqy[2:(maxFreqOffsetIdx+1)]), 0, freqy[2:(maxFreqOffsetIdx+1)])
  } else {
    oFreq <- 0
  }

  # determine block lengths
  blockx <- blockStartIdx(n = nx, blockSize = blockSizex, overlap = overlap)
  blocky <- blockStartIdx(n = ny, blockSize = blockSizey, overlap = overlap)

  # calculate slepians to (hopefully?) speed up spec.mtm a smidge.
  slepy <- dpss(n = blockSizey, k = k, nw = nw, returnEigenvalues = FALSE)$v
  slepx <- dpss(n = blockSizex, k = k, nw = nw, returnEigenvalues = FALSE)$v

  # initial storage based on freqRangeIdx -/+ maxFreqOffsetIdx
  yky <- array(NA_complex_, dim = c(nFreqRangeIdx, k, blocky$nBlock))
  ykx <- array(NA_complex_, dim = c(nFreqRangeIdx + nOffsets - 1, k, blockx$nBlock, dim(x)[2]))

  ### put this loopiness into another function? Check copying in functions... when will
  # R make copies?
  for (i in 1:blocky$nBlock){
    for (j in 1:ncol(x)){
      tmp <- eigenCoef(x[ blockx$startIdx[i]:(i*blockSizex), j]
                                 , nw = nw, k = k, nFFT = nFFTx, centre = centre
                                 , deltat = dtx, dpssIN = slepx)
      # add negative frequencies if needed:
      if (freqRangeIdx[1] - maxFreqOffsetIdx <= 0){
        ykx[, , i, j] <- eigenCoefWithNegYk(tmp, freqRangeIdx, maxFreqOffsetIdx)
      } else {
        ykx[, , i, j] <- tmp[(freqRangeIdx[1]-maxFreqOffsetIdx):(freqRangeIdx[2]-maxFreqOffsetIdx)
                             , ]
      }
    }
    yky[, , i] <- eigenCoef(y[ blocky$startIdx[i]:(i*blockSizey), 1]
                            , nw = nw, k = k, nFFT = nFFTy, centre = centre
                            , deltat = dty, dpssIN = slepy)[freqRangeIdx[1]:freqRangeIdx[2], ]
  }

  if (is.null(decorrelate)){
    coh <- array(0, dim = c(nOffsets, nFreqRangeIdx, ncol(x)))
    for (i in 1:blocky$nBlock){
      for (j in 1:ncol(x)){
        coh[, , j] <- coh[, , j] + mscFromEigenHelper(yk1 = yky[, , i]
                                                      , yk2 = ykx[, , i, j]
                                                      , k = k
                                                      , nOffsets = nOffsets)
      }
    }
    coh <- coh / blocky$nBlock

    # manipulate the coherence matrix to enforce:
    # 1) zeroFreqOffset
    # 2) deadZone (from filtering)
    if (forceZeroOffset){
      coh[maxFreqOffsetIdx+1, , ] <- 999
    }

    ##### TODO: add checks - need freqRange to start at 0 I believe for this to work the way you would expect ...
    ### Fix this function to work in all cases regardless, needs some messing around
    if (!is.null(deadBand)){
      for (j in 1:ncol(x)){
        coh[, , j] <- fixDeadBand(msc = coh[, , j]
                           , zeroOffsetIdx = maxFreqOffsetIdx + 1
                           , freqRangeIdx = freqRangeIdx, band = deadBand[2], df = df
                           , replaceWith = -999)

      }
    }

    mscLev <- qnorm(cohSigLev, mean = 0, sd = 1 / sqrt(blocky$nBlock))

    coh.ind <- array(NA_integer_, dim = dim(coh))

    if (maxFreqOffsetIdx == 0){
      coh.ind[, , ] <- 1
    } else {
      # think about how to optimize this - specifically, memory usage is likely inefficient
      for(j in 1:ncol(x)){
        coh.ind[, , j] <- matrix(.Fortran("msc2indicator"
                                          , msc = as.double(coh[, , j])
                                          , nrow = as.integer(dim(coh)[1])
                                          , ncol = as.integer(dim(coh)[2])
                                          , ind = integer(dim(coh)[1]*dim(coh)[2])
                                          , level = as.double(mscLev)
                                          , nOff = as.integer(nOffAllowed))$ind
                                 , nrow = dim(coh)[1], ncol = dim(coh)[2])
      }
    }

    coh.ind[coh.ind > 0] <- 1

    # out <- .Fortran("tf"
    #                 , H
    #                 , yk1
    #                 , yk2
    #                 , cohInd
    #                 , nrow1
    #                 , nrow2
    #                 , npred
    #                 , nBlocks
    #                 , nOffsets
    #                 , nUniqOffsets
    #                 , k)
  } else {

  }



  info <- list(namey = namey, namex = namex
               , method = method, nw = nw, k = k, , nFFTy = nFFTy, nFFTx = nFFTx
               , centre = centre, dty = dty, dtx = dtx, dtRatio = dtRatio
               , blockSizey = blockSizey, blockSizex = blockSizex
               , overlap = overlap
               , nBlocks = blocky$nBlock
               , freqRange = freqRange, freqRangeIdx = freqRangeIdx
               , nFreqRangeIdx = nFreqRangeIdx
               , maxFreqOffset = maxFreqOffset, maxFreqOffsetIdx = maxFreqOffsetIdx
               , convertMscToNormal = TRUE)
}

tfEigen <- function(yk1, yk2, cohInd){

}
