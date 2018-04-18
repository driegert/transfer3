# coherence based methods:

#' @details Only calculates the forward coherence currently.
#' @export
coherence <- function(x, y, method = c("aveOfMsc", "mscOfAve", "mscOfAveCoh")
                      , nw = 4, k = 7, nFFTx = NULL
                      , dtx = 1, dty = 1, blockSizex = NULL, overlap = 0
                      , freqRange = NULL, maxFreqOffset = 0
                      , convertMscToNormal = FALSE
                      , namex = "x", namey = "y"){

  dtRatio <- dty / dtx
  nx <- length(x)
  ny <- length(y)

  if (is.null(blockSizex)){
    blockSizex <- nx
  }

  # get correct blocksizes:
  stopifnot(dtRatio >= 1 # series y should be sampled as or more quickly than series x
            , blockSizex <= length(x) # blockSize can't be larger than the length of the series
            , dtRatio - trunc(dtRation) == 0) # need the ratio to be an integer for this to work properly

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
    freqRange <- c(0, 1/(2*dtx))
    freqRangeIdx <- c(1, nfreqx)
  }

  maxFreqOffsetIdx <- ceiling(maxFreqOffset / df)

  blockStartx <- blockStartIdx(n = nx, blockSize = blockSizex, overlap = overlap)
  ### start here:
  blockStarty <- list(startIdx = (blockStartx$increment - 1) * dtRatio + 1
                      , increment = blockStartx$increment
                      , nBlock = blockStartx$nBlock)

}

# the assumption here is that the df's are the same in yk1 and yk2, which means nFFT1 = a*nFFT2 (a = some integer)
aveOfMsc <- function(yk1, yk2, freqRangeIdx, maxOffsetIdx){
  if (freqRangeIdx[2] + maxOffsetIdx > nrow(yk2)){
    stop("Not enough indices to accommodate the highest offset from the largest centre frequency.")
  }

  nfreq1 <- freqRangeIdx[1] - freqRangeIdx[2] + 1
  nOffsets <- 2*maxOffsetIdx + 1

  if (freqRangeIdx[1] - maxOffsetIdx <= 0){
    yk2.f <- rbind(Conj(yk2[rev(2:(maxOffsetIdx - freqRangeIdx[1] + 2)), ])
                   , yk2[1:(freqRangeIdx[2] + maxOffsetIdx), ])
    out <- .Fortran("cohMsc"
                    , coh = complex(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2] ])
                    , yk2 = as.complex(yk2.f)
                    , nfreq1 = as.integer(nfreq1)
                    , nfreq2 = as.integer(nrow(yk2.f))
                    , k = as.integer(ncol(yk1))
                    , nOffsets = as.integer(nOffsets))
  } else {
    out <- .Fortran("cohMsc"
                    , coh = complex(nfreq1 * nOffsets)
                    , yk1 = as.complex(yk1[ freqRangeIdx[1]:freqRangeIdx[2] ])
                    , yk2 = as.complex(yk2[ freqRangeIdx[1]:freqRangeIdx[2] ])
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
