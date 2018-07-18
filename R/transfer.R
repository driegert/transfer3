# transfer function related code:

#' Calculates some transfer functions.
#'
#' @param y A single column \code{data.frame} containing the response time-series.
#' @param x A \code{data.frame} containing the predictor time-series.
#' @param decorrelate A vector containing the names of the predictors (columns of x)
#' in the order of the iterative regression.
#' @param nOffAllowed An \code{integer} representing the number of offsets, \emph{including
#' the zero-offset}, to be used in the transfer function estimation.
#' @param deadband A vector of two values indicating the start and end of the deadband to remove from the coherence
#' matrix - this doesn't work if freqRange[1] != 0 - also, you're going to get errors if you run this when
#' maxFreqOffset == 0.
#' @export
tf <- function(y, x, nw = 4, k = 7, nFFTy = NULL, centre = "none"
               , dty = 1, dtx = 1, blockSizey = NULL, overlap = 0
               , cohSigLev = 0.99, nOffAllowed = 1
               , forceZeroOffset = TRUE
               , freqRange = NULL, maxFreqOffset = 0, deadBand = NULL
               , decorrelate = NULL, nDecorOffAllowed = NULL){
  ## check ALL the parameters here...
  # 1) make sure decorrelate contains column names of x
  # 2) dty > dtx
  # 3) deadBand has c(0, <something>)
  # 4) if deadBand != NULL, freqRange = c(0, <something>)
  # add some checks for NA's here
  #################
  #
  #
  #
  #################

  ## TODO: need to double check the number of frequencies being used:
  ## especially wrt to decorrelation.  nDecorFreqIdx might not be right if we don't start at f = 0 ... due to
  # the negative frequencies...

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
    nFFTy <- 2^ceiling(log2(ny) + 1)
  }
  nFFTx <- dtRatio * nFFTy

  freqy <- posFreq(nFFTy, dty)
  freqx <- posFreq(nFFTx, dtx)

  # frequency spacing
  df <- 1/(nFFTy * dty)

  # find the number of indices in the band [0, W]
  w <- nw / (dty * blockSizey)
  nwIdx <- floor(w / df)

  nfreqy <- length(freqy)
  nfreqx <- length(freqx)

  if (is.null(freqRange)){
    freqRange <- c(0, 1/(2*dty))
    freqRangeIdx <- c(1, nfreqy)
  } else {
    freqRangeIdx <- c(max(1, floor(freqRange[1] / df)), min(floor(freqRange[2] / df), nfreqy))
  }

  nFreqRangeIdx <- freqRangeIdx[2] - freqRangeIdx[1] + 1

  maxFreqOffsetIdx <- ceiling(maxFreqOffset / df)

  # Different frequency ranges based on if we're decorrelating or not.
  # keeps more rows in the eigencoefficient calculations if decorrelating
  if (!is.null(decorrelate)){
    decorFreqRangeIdx <- freqRangeIdx + c(-maxFreqOffsetIdx, maxFreqOffsetIdx)
    decorFreqRangeIdx[1] <- max(1, decorFreqRangeIdx[1])
  } else {
    decorFreqRangeIdx <- freqRangeIdx
  }

  nOffsets <- 2*maxFreqOffsetIdx + 1
  nDecorFreqRangeIdx <- decorFreqRangeIdx[2] - decorFreqRangeIdx[1] + 1 ### TODO: deal with this...

  if (!is.null(decorrelate) && is.null(nDecorOffAllowed)){
    nDecorOffAllowed <- nOffAllowed
  }

  if (freqRangeIdx[1] != 1 & !is.null(decorrelate)){
    stop("I don't think this works - check the TODO comments in tf() and remove this stop() line if it seems okay")
  }

  if (nOffsets > 1){
    oFreq <- c(rev(-1*freqy[2:(maxFreqOffsetIdx+1)]), 0, freqy[2:(maxFreqOffsetIdx+1)])
  } else {
    oFreq <- 0
  }

  # determine block lengths
  blockx <- blockStartIdx(n = nx, blockSize = blockSizex, overlap = overlap)
  blocky <- blockStartIdx(n = ny, blockSize = blockSizey, overlap = overlap)

  # calculate slepians to (hopefully?) speed up spec.mtm a smidge.
  slepy <- multitaper::dpss(n = blockSizey, k = k, nw = nw, returnEigenvalues = FALSE)$v
  slepx <- multitaper::dpss(n = blockSizex, k = k, nw = nw, returnEigenvalues = FALSE)$v

  # initial storage based on freqRangeIdx -/+ maxFreqOffsetIdx
  yky <- array(NA_complex_, dim = c(nFreqRangeIdx, k, blocky$nBlock))
  ykx <- array(NA_complex_, dim = c(nDecorFreqRangeIdx + nOffsets - 1, k, blockx$nBlock, dim(x)[2]))

  ### put this loop into another function? Check copying in functions... when will
  # R make copies?
  for (i in 1:blocky$nBlock){
    for (j in 1:ncol(x)){
      tmp <- eigenCoef(x[blockx$startIdx[i]:(i*blockSizex), j]
                                 , nw = nw, k = k, nFFT = nFFTx, centre = centre
                                 , deltat = dtx, dpssIN = slepx)
      # add negative frequencies if needed:
      if (decorFreqRangeIdx[1] - maxFreqOffsetIdx <= 0){
        ykx[, , i, j] <- eigenCoefWithNegYk(tmp, decorFreqRangeIdx, maxFreqOffsetIdx)
      } else {
        ykx[, , i, j] <- tmp[(decorFreqRangeIdx[1]-maxFreqOffsetIdx):(decorFreqRangeIdx[2]+maxFreqOffsetIdx)
                             , ]
      }
    }
    yky[, , i] <- eigenCoef(y[ blocky$startIdx[i]:(i*blockSizey), 1]
                            , nw = nw, k = k, nFFT = nFFTy, centre = centre
                            , deltat = dty, dpssIN = slepy)[freqRangeIdx[1]:freqRangeIdx[2], ]
  }

  ## Decorrelate the inputs before calculating the transfer functions
  # store eigencoefficients and return? No - return decorrelated time-series.
  # no - return the eigencoefficients.  We can get the time-series from those by
  # using complex demodulates and narrow-band frequency information if needed?
  # I don't *think* we actually need the time-series, in HC, we use the impulse responses
  # as a representation of risk..

  # indices for decorrelation are becoming a pain in my ass.  Book keeping is the absolute *worst* for this problem.

  if (!is.null(decorrelate)){
    #### TODO: deal with this...
    # decorCentFreqIdx <- (dim(ykx)[1] - nDecorFreqRangeIdx + 1):(dim(ykx)[1] - maxFreqOffsetIdx) # old, and wrong?
    decorCentFreqIdx <- (maxFreqOffsetIdx - freqRangeIdx[1] + 2):(dim(ykx)[1] - maxFreqOffsetIdx)

    ####################
    # Calculate the L's
    ####################
    ####################
    nameOrder <- rep(-1, length(decorrelate))
    for (i in 1:length(decorrelate)){
      nameOrder[i] <- which(decorrelate[i] == names(x))
    }

    for (l in 2:length(decorrelate)){
      coh <- array(0, dim = c(nOffsets, nDecorFreqRangeIdx, l-1))

      for (i in 1:blockx$nBlock){
        for (j in 1:(l-1)){
          coh[, , j] <- coh[, , j] + mscFromEigenHelper(yk1 = ykx[decorCentFreqIdx, , i, nameOrder[l]]
                                                        , yk2 = ykx[, , i, nameOrder[l-j]]
                                                        , k = k
                                                        , nOffsets = nOffsets)
        }
      }
      coh <- coh / blockx$nBlock

      # deal with potential infinities:
      coh[is.infinite(coh) & coh < 0] <- min(coh[!is.infinite(coh)])
      coh[is.infinite(coh) & coh > 0] <- max(coh[!is.infinite(coh)])

      # manipulate the coherence matrix to enforce:
      # 1) zeroFreqOffset
      # 2) deadZone (from filtering)

      #1)
      if (forceZeroOffset){
        ### START HERE - MAY 13, 2018
        centreBandW <- (maxFreqOffsetIdx+1-floor(nwIdx / 2)):(maxFreqOffsetIdx+1+floor(nwIdx/2))
        if (any(centreBandW < 1)){
          coh[, , ] <- -999
        } else {
          coh[centreBandW, , ] <- -999
        }
        coh[maxFreqOffsetIdx+1, , ] <- 999
      }

      #2)
      ##### TODO: add checks - need freqRange and deadBand to start at 0 I believe for
      ### this to work the way you would expect ...
      ### Fix this function to work in all cases regardless, needs some messing around
      ### Might just want to write this in Fortran ...
      if (!is.null(deadBand)){
        for (j in 1:(l-1)){
          coh[, , j] <- fixDeadBand(msc = coh[, , j]
                                    , zeroOffsetIdx = maxFreqOffsetIdx + 1
                                    , freqRangeIdx = freqRangeIdx, band = deadBand[2], df = df
                                    , replaceWith = -999)
        }
      }

      mscLev <- qnorm(cohSigLev, mean = 0, sd = 1 / sqrt(blockx$nBlock))

      coh.ind <- array(NA_integer_, dim = dim(coh))

      if (maxFreqOffsetIdx == 0){
        coh.ind[, , ] <- 1
      } else {
        # think about how to optimize this - specifically, memory usage is likely inefficient
        # "modify in place" would be great?
        for(j in 1:(l-1)){
          coh.ind[, , j] <- matrix(.Fortran("msc2indicator"
                                            , msc = as.double(coh[, , j])
                                            , nrow = as.integer(dim(coh)[1])
                                            , ncol = as.integer(dim(coh)[2])
                                            , ind = integer(dim(coh)[1]*dim(coh)[2])
                                            , level = as.double(mscLev)
                                            , nOff = as.integer(nDecorOffAllowed))$ind
                                   , nrow = dim(coh)[1], ncol = dim(coh)[2])
        }
      }

      coh.ind[coh.ind > 0] <- 1

      # TODO: fix this for speed - Fortran maybe?
      HcolIdx <- which(apply(coh.ind, c(1,3), sum) > 0, arr.ind = TRUE) #number of 1's in each row, each predictor
      # "col" indicates the predictor
      nUniqueOffsets <- dim(HcolIdx)[1]

      H <- matrix(.Fortran("tf"
                           , H = complex(nUniqueOffsets * nDecorFreqRangeIdx)
                           , yk1 = as.complex(ykx[decorCentFreqIdx, , , nameOrder[l]])
                           , yk2 = as.complex(ykx[, , , nameOrder[1:(l-1)]])
                           , cohInd = as.integer(coh.ind)
                           , nrow1 = as.integer(nDecorFreqRangeIdx)
                           , nrow2 = as.integer(nrow(ykx))
                           , npred = as.integer(l-1)
                           , nBlocks = as.integer(blockx$nBlock)
                           , nOffsets = as.integer(nOffsets)
                           , nUniqOffsets = as.integer(nUniqueOffsets)
                           , k = as.integer(k))$H
                  , nrow = nDecorFreqRangeIdx, ncol = nUniqueOffsets)

      Hinfo <- HcolInfo(maxFreqOffsetIdx = maxFreqOffsetIdx
                        , df = df, npred = l-1
                        , HcolIdx = HcolIdx
                        , predNames = names(x)[nameOrder[1:(l-1)]])


      xLfitted <- eigenCoefFit(H, Hinfo, ykx[, , , nameOrder[1:(l-1)], drop = FALSE]
                               , names(x)[nameOrder[1:(l-1)]], decorCentFreqIdx)

      # appends the negative frequencies on to the top if needed to line up with ykx properly.
      if (decorFreqRangeIdx[1] - maxFreqOffsetIdx <= 0){
        xLfitted <- eigenCoefWithNegYk(yk2 = xLfitted
                                       # , freqRangeIdx = c(1, dim(xLfitted)[1])
                                       , freqRangeIdx = freqRangeIdx
                                       , maxFreqOffsetIdx = maxFreqOffsetIdx, multPred = TRUE)
      }

      ykx[1:dim(xLfitted)[1], , , nameOrder[l]] <- ykx[1:dim(xLfitted)[1], , , nameOrder[l], drop = FALSE] - xLfitted
    }

    # done right up there ^^ just a couple lines up ^^
    # # Replace the negative frequencies if needed: ## Again - this only works if freqRangeIdx[1] == 1 I think...
    # # also trims off the higher frequencies that aren't needed for the last transfer function estimation.
    # if (decorFreqRangeIdx[1] - maxFreqOffsetIdx <= 0){
    #   ykx <- eigenCoefWithNegYk(ykx[decorCentFreqIdx, , , , drop = FALSE], freqRangeIdx
    #                             , maxFreqOffsetIdx, multPred = TRUE)
    # } else {
    #   ykx <- ykx[decorCentFreqIdx[1:nFreqRangeIdx], , , , drop = FALSE]
    # }
  }

  ## It's business time.  I know it's business time because it's time to start
  # the calculating of the H's
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

  # deal with potential infinities:
  coh[is.infinite(coh) & coh < 0] <- min(coh[!is.infinite(coh)])
  coh[is.infinite(coh) & coh > 0] <- max(coh[!is.infinite(coh)])

  # manipulate the coherence matrix to enforce:
  # 1) zeroFreqOffset
  # 2) deadZone (from filtering)

  #1)
  if (forceZeroOffset){
    centreBandW <- (maxFreqOffsetIdx+1-floor(nwIdx / 2)):(maxFreqOffsetIdx+1+floor(nwIdx/2))
    if (any(centreBandW < 1)){
      coh[, , ] <- -999
    } else {
      coh[centreBandW, , ] <- -999
    }
    coh[maxFreqOffsetIdx+1, , ] <- 999
  }

  #2)
  ##### TODO: add checks - need freqRange and deadBand to start at 0 I believe for
  ### this to work the way you would expect ...
  ### Fix this function to work in all cases regardless, needs some messing around
  ### Might just want to write this in Fortran ...
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

  # browser()
  #### ADDED && forceZeroOffset on July 18, 2018
  if (maxFreqOffsetIdx == 0 && forceZeroOffset){
    coh.ind[, , ] <- 1
    # ADDED this entire else if() statement on July 18, 2018 as well.
  } else if (maxFreqOffsetIdx == 0 && !forceZeroOffset){
    coh.ind[, , ] <- 0
    coh.ind[coh > mscLev] <- 1
  } else {
    # think about how to optimize this - specifically, memory usage is likely inefficient
    # "modify in place" would be great?
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

  # TODO: fix this for speed - Fortran maybe?
  HcolIdx <- which(apply(coh.ind, c(1,3), sum) > 0, arr.ind = TRUE)
  nUniqueOffsets <- dim(HcolIdx)[1]

  H <- matrix(.Fortran("tf"
                       , H = complex(nUniqueOffsets * nFreqRangeIdx)
                       , yk1 = as.complex(yky)
                       , yk2 = as.complex(ykx)
                       , cohInd = as.integer(coh.ind)
                       , nrow1 = as.integer(nrow(yky))
                       , nrow2 = as.integer(nrow(ykx))
                       , npred = as.integer(ncol(x))
                       , nBlocks = as.integer(blockx$nBlock)
                       , nOffsets = as.integer(nOffsets)
                       , nUniqOffsets = as.integer(nUniqueOffsets)
                       , k = as.integer(k))$H
              , nrow = nFreqRangeIdx, ncol = nUniqueOffsets)

  Hinfo <- HcolInfo(maxFreqOffsetIdx = maxFreqOffsetIdx
                    , df = df, npred = dim(x)[2]
                    , HcolIdx = HcolIdx
                    , predNames = names(x))

  info <- list(namey = names(y), namex = names(x)
               , nw = nw, k = k, nFFTy = nFFTy, nFFTx = nFFTx
               , centre = centre, dty = dty, dtx = dtx, dtRatio = dtRatio
               , blockSizey = blockSizey, blockSizex = blockSizex
               , overlap = overlap
               , nBlocks = blocky$nBlock
               , freqRange = freqRange, freqRangeIdx = freqRangeIdx
               , nFreqRangeIdx = nFreqRangeIdx
               , maxFreqOffset = maxFreqOffset, maxFreqOffsetIdx = maxFreqOffsetIdx
               , df = df, convertMscToNormal = TRUE)

  if (is.null(decorrelate)){
    list(H = H, Hinfo = Hinfo, info = info, ykx = NULL)
  } else {
    list(H = H, Hinfo = Hinfo, info = info
         , ykx = ykx[decorCentFreqIdx[1]:dim(ykx)[1], , , , drop = FALSE])
  }
}


#' Impulse response from transfer functions
#'
#' Inverse Fourier transforms the transfer function to give the impulse response
#'
#' @param H A \code{matrix} where each column contains a transfer function - tf()$H would be appropriate.
#' @param n An \code{integer} indicating the half-length of the impulse response to return, L = 2*n+1.
#' Normally this would be the \code{blockSize2} agrument used in tf().
#' @param realPart A \code{logical} indicating whether to return only the real part of the impulse response
#' (default = TRUE) or return the complex-valued impulse response (FALSE)
#'
#' @export
ir <- function(H, n, realPart = TRUE){
  stopifnot(is.matrix((H)) | n > ncol(H))
  # "negative" frequencies in the top half of the array
  # ^^ these are conjugated first and reversed (should be conjugate symmetric after)
  Hfull <- rbind(H, Conj(H[(nrow(H) - 1):2, , drop = FALSE]))

  # should be real-valued with real-valued data as inputs.
  if (realPart){
    h <- Re(mvfft(Hfull, inverse = TRUE) / nrow(Hfull))
  } else {
    h <- mvfft(Hfull, inverse = TRUE) / nrow(Hfull)
  }

  n <- n+1
  # put the impulse response in the correct place in the array in order to use in a convolution - i.e., filter()
  ind <- c((nrow(h)-n+2):nrow(h), 1:n)

  list(h = h[ind, , drop = FALSE], n = n-1, realPart = realPart)
}


#' Spectrum prediction
#' Predicts the spectrum based on an estimated transfer function and new data provided.
#'
#' @param H A \code{list} returned by tf() containing H, Hinfo, and info.
#' @param d2 A \code{data.frame} containing the new data.  Must have the same column names as the original
#' data used in the transfer function estimation.
#'
#' @export
specPredict <- function(H, d2){
  predNames <- names(d2)
  spec <- list()
  for (i in 1:length(predNames)){
    spec[[predNames[i]]] <- multitaper::spec.mtm(d2[, i], nw = H$info$nw, k = H$info$k
                                                 , deltat = H$info$dtx
                                                 , nFFT = H$info$nFFTx
                                                 , center = 'none'
                                                 , returnInternals = TRUE, plot = FALSE)$spec
  }

  fullFreqRange <- H$info$freqRangeIdx[1]:H$info$freqRangeIdx[2]

  sRecon <- rep(0, length(fullFreqRange))

  for (i in 1:dim(H$Hinfo)[1]){
    offIdx <- fullFreqRange + H$Hinfo$idxOffset[i]
    offIdx[offIdx <= 0] <- abs(offIdx[offIdx <= 0]) + 2
    sRecon <- sRecon + (abs(H$H[, i])^2) * spec[[ H$Hinfo$predictor[i] ]][offIdx]
  }

  sRecon
}

#' Predict the Time Series from the predictors and impulse response
#'
#' @param newdata A \code{data.frame} containing the predictors.  Column names must match what is
#' present in Hinfo.
#' @param ir A \code{matrix} containing the impulse responses corresponding to the rows in Hinfo.
#' @param Hinfo A piece of what is returned by tf().
#' @param info Another piece of what is returned by tf().
#'
#' @export
predictTs <- function(newdata, ir, Hinfo, info){
  N <- nrow(newdata)
  yhat <- rep(0, N)
  # nFlt <- nrow(ir)
  nTrim <- ir$n

  # # standardize the input if necessary.
  # if (info$standardize){
  #   newdata <- df.std(newdata)
  # }

  for (i in 1:nrow(Hinfo)){
    ## Should this be 2*pi*(1:N) or 2*pi*(0:(N-1)) ?
    ## is the negative necessary? I guess it doesn't matter since cos(x) is an even function...
    # 0:(N-1) inside cos() due to fft ranging from 0:(N-1) as well.
    if (Hinfo$idxOffset[i] == 0){
      tmp <- filter.twoSided(newdata[, Hinfo$predictor[i]], ir$h[, Hinfo$Hcolumn[i]])
    } else {
      tmp <- filter.twoSided(2*cos(-2*pi*(0:(N-1))*(Hinfo$idxOffset[i]/info$nFFTy)) * newdata[, Hinfo$predictor[i]]
                     , ir$h[, Hinfo$Hcolumn[i]])  ####### CHECK THIS - cos() - I'm unsure.
    }

    yhat <- yhat + tmp # this adds NA's at the end, also has NA's at the beginning - I need to look at
    # zFilter again for how this works by convoling using FFT's ...
  }

  yhat
}


#' @export
HcolInfo <- function(maxFreqOffsetIdx, df, npred, HcolIdx, predNames){
  ### Work from here:
  # book keeping - which frequencies go where?
  ####
  # copied directly from transfer2::tf()
  # - April 27, 2018
  ###
  # we need to obtain the value of the offsets to be used
  # - probably actualy index and also in terms of frequency
  offIdxFromCent <- (-maxFreqOffsetIdx):maxFreqOffsetIdx
  offFreqFromCent <- offIdxFromCent * df

  # create a data.frame of info for the transfer function matrix.
  # can probably do this in one line...
  for (i in 1:npred){
    if (i == 1){
      curVal <- 1:length(HcolIdx[ HcolIdx[, "col"] == i, "row" ])
      Hinfo <- data.frame(predictor = as.character(predNames[i])
                          , Hcolumn = curVal
                          , freqOffset = offFreqFromCent[ HcolIdx[ HcolIdx[, "col"] == i, "row" ] ]
                          , idxOffset = offIdxFromCent[ HcolIdx[ HcolIdx[, "col"] == i, "row" ] ]
                          , stringsAsFactors = FALSE)
    } else {
      curVal <- tail(curVal, 1) + 1:length(HcolIdx[ HcolIdx[, "col"] == i, "row" ])
      Hinfo <- rbind(Hinfo,data.frame(predictor = as.character(predNames[i])
                                      , Hcolumn = curVal
                                      , freqOffset = offFreqFromCent[ HcolIdx[ HcolIdx[, "col"] == i, "row" ] ]
                                      , idxOffset = offIdxFromCent[ HcolIdx[ HcolIdx[, "col"] == i, "row" ] ]
                                      , stringsAsFactors = FALSE)
      )
    }
    # hColNames <- c(hColNames
    #                , paste0(names(d2)[i], "..", offFreqFromCent[hIdx[ hPredBreak[i]:(hPredBreak[i+1] - 1) ]]))
  }

  Hinfo
}
