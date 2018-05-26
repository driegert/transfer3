# testing the Fortran code to make sure it makes sense!

#' @useDynLib transfer3
test.cohMsc <- function(){
  set.seed(19840331)
  yk1 <- matrix(complex(real = rnorm(12), imaginary = rnorm(12)), ncol = 3, nrow = 4)
  yk2 <- matrix(complex(real = rnorm(24), imaginary = rnorm(24)), ncol = 3, nrow = 8)

  nOffsets <- 5 # 2*2 + 1
  nfreq1 <- 4
  nfreq2 <- 8
  k <- 3

  out <- .Fortran("cohMsc"
                  , coh = double(nOffsets * nfreq1)
                  , yk1 = as.complex(yk1)
                  , yk2 = as.complex(yk2)
                  , nfreq1 = as.integer(nfreq1)
                  , nfreq2 = as.integer(nfreq2)
                  , k = as.integer(k)
                  , nOffsets = as.integer(nOffsets))

  matrix(out$coh, nrow = nOffsets, ncol = nfreq1)
}

# matrices in and out as you would expect.
#' @useDynLib transfer3
test.testMatrix <- function(){
  x <- matrix(1:12, nrow = 4)
  out <- .Fortran("testMatrix", x = as.integer(x))

  y <- matrix(out$x, nrow = 4)

  y
}

verify.cohMsc <- function(){
  set.seed(19840331)
  yk1 <- matrix(rnorm(12), ncol = 3, nrow = 4)
  yk2 <- matrix(rnorm(24), ncol = 3, nrow = 8)

  nOffsets <- 5 # 2*2 + 1
  nfreq1 <- 4
  nfreq2 <- 8
  k <- 3

  s1 <- apply(abs(yk1)^2, 1, sum)
  s2 <- apply(abs(yk2)^2, 1, sum)

  coh <- matrix(NA, nrow = nOffsets, ncol = nfreq1)

  for (i in 1:nOffsets){
    coh[i, ] <- abs(apply(yk1 * yk2[i:(nfreq1+i-1), ], 1, sum))^2 / (s1 * s2[i:(nfreq1+i-1)])
  }

  coh
}


test.r.cohMsc <- function(){
  x <- rnorm(100)
  y <- rnorm(100)

  nFFT <- 256

  coh2 <- transfer2::coherence(d1 = x, d2 = y, ndata = 100, ndata2 = 100
                               , blockSize = 100, blockSize2 = 100
                               , nFFT = 256, nFFT2 = nFFT, dt = 1, dt2 = 1)

}

test.designMatSetup <- function(){
  k <- 3
  npred <- 1
  nOff <- 2
  nOffsets <- 7 # maxFreqOffsetIdx = 3 in this case
  yk1 <- matrix(complex(real = rnorm(6*k), imaginary = rnorm(6*k)), ncol = k)
  nrow2 <- 12
  yk2 <- array(complex(real = rnorm(12*k*npred), imaginary = rnorm(12*k*npred)), dim = c(nrow2, k, npred))

  coh <- matrix(0, ncol = nrow(yk1), nrow = nOffsets)
  coh[4, ] <- 1
  coh[2, 1] <- 1
  coh[5, 2] <- 1
  coh[2, 5] <- 1
  coh[7, 6] <- 1

  nUniqOff <- sum(apply(coh, 1, sum) > 0)

  out <- .Fortran("tf", H = complex(nrow(yk1)*nUniqOff)
                  , yk1 = as.complex(yk1)
                  , yk2 = as.complex(yk2)
                  , cohInd = as.integer(coh)
                  , nrow1 = as.integer(nrow(yk1))
                  , nrow2 = as.integer(nrow(yk2))
                  , npred = as.integer(npred)
                  , nBlocks = as.integer(1)
                  , nOffsets = as.integer(nOffsets)
                  , nUniqOffsets = as.integer(nUniqOff)
                  , k = as.integer(k))
}

test.tf <- function(){
  print("No offsets case.  Exponential transfer function.")
  x <- arima.sim(n = 100, model = list(ar = c(0.9, -0.2)))
  X <- eigenCoef(x)

  freq <- seq(0, 0.5, length.out = 129)
  H <- (freq+1)^4

  Y <- X
  for (i in 1:ncol(X)){
    Y[, i] <- H * X[, i] + complex(real = rnorm(dim(X)[1]), imaginary = rnorm(dim(X)[1]))
  }

  nUniqOff <- 1
  out <- .Fortran("tf", H = complex(nrow(Y)*nUniqOff)
                  , yk1 = as.complex(Y)
                  , yk2 = as.complex(X)
                  , cohInd = as.integer(matrix(1, nrow = 1, ncol = 129))
                  , nrow1 = as.integer(nrow(Y))
                  , nrow2 = as.integer(nrow(X))
                  , npred = as.integer(1)
                  , nBlocks = as.integer(1)
                  , nOffsets = as.integer(1)
                  , nUniqOffsets = as.integer(nUniqOff)
                  , k = as.integer(ncol(X)))

  plot(Re(out$H), type='l', ylim = c(min(Im(out$H)), max(Re(out$H))))
  lines(H, col = 'black', lty = 2)
  lines(Im(out$H), col = 'blue')
  lines(Im(H), col = 'blue', lty = 2)
}


test.r.tf <- function(){
  x <- data.frame(x1 = arima.sim(n = 200, model = list(ar = c(0.9, -0.2)))
                  , x2 = arima.sim(n = 200, model = list(ar = c(0.7, -0.2)))
                  , x3 = arima.sim(n = 200, model = list(ar = c(0.5, -0.1))))
  y <- data.frame(y = arima.sim(n = 100, model = list(ar = c(0.9, -0.2))))

  nw = 4; k = 7; nFFTy = NULL; centre = "none"
  dty = 2; dtx = 1; blockSizey = NULL; overlap = 0
  cohSigLev = 0.9; nOffAllowed = 4 ### careful here!
  forceZeroOffset = TRUE
  freqRange = NULL; maxFreqOffset = .02
  decorrelate = NULL
  deadBand = c(0, 0.004)
  forceZeroOffset = TRUE
  maxFreqOffset = 0.1
  decorrelate = c("x3", "x1", "x2")


  tfEst <- tf(y = y, x = x, nw = nw, k = k, dty = dty, dtx = dtx, cohSigLev = 0.9
              , nOffAllowed = nOffAllowed
              , forceZeroOffset = TRUE, maxFreqOffset = 0.02, decorrelate = decorrelate
              , deadBand = c(0, 0.004))


  # check if the decorrelation actually even did anything...
  x.spec <- lapply(x, eigenCoef, nw = nw, k = k, nFFT = "default"
                   , deltat = dtx)
  diff1 <- tfEst$ykx[1:140, , , 1] - x.spec$x1[1:140, ]
  for (i in 1:k){ plot(abs(diff1[, i]), type='l') }

  specRecon <- specPredict(tfEst, d2 = x)

  par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
  sy <- spec.mtm(y, deltat = 2, dtUnits = "second", nw = nw, k = k, main = "")
  lines(sy$freq, specRecon, col = 'blue')

  # X <- eigenCoef(x)
  #
  # freq <- seq(0, 0.5, length.out = 129)
  # H <- (freq+1)^4
  #
  # Y <- X
  # for (i in 1:ncol(X)){
  #   Y[, i] <- H * X[, i] + complex(real = rnorm(dim(X)[1]), imaginary = rnorm(dim(X)[1]))
  # }

  ir <- ir(H = tfEst$H, n = 20)

}

test.twosided <- function(){
  filter <- c(1, 2, 3, 4, 5)

  series <- 1:25

  out <- .Fortran("twoSidedFilter"
                  , filter = as.double(filter)
                  , series = as.double(series)
                  , seriesOut = double(length(series))
                  , hTrim = as.integer(floor(length(filter) / 2))
                  , n = as.integer(length(series)))
  out$seriesOut[which(out$seriesOut == -999.99)] <- NA

  rout <- filter(series, filter = filter, sides = 2)
}


test.decor.coh <- function(){
  png("~/school_lab/bin/testing/transfer3/plots/beforeDecorrelate.png"
      , height = 6, width = 8, res = 200, units = "in")
  par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
  fields::image.plot(x = (1e3)*df*(0:(dim(coh1)[2]-1)), y = (1e3)*df*((-maxFreqOffsetIdx):maxFreqOffsetIdx)
                     , z=t(coh1[,,1])
                     , zlim = range(c(coh1, coh2))
                     , xlab = "Central Frequency (mHz)"
                     , ylab = "Offset Frequency (mHz)"
                     , legend.lab = expression(paste("Coherence (", sigma, " from mean)")))
  dev.off()

  png("~/school_lab/bin/testing/transfer3/plots/afterDecorrelate.png"
      , height = 6, width = 8, res = 200, units = "in")
  par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
  fields::image.plot(x = (1e3)*df*(0:(dim(coh2)[2]-1)), y = (1e3)*df*((-maxFreqOffsetIdx):maxFreqOffsetIdx)
                     , z=t(coh2[,,1])
                     , zlim = range(c(coh1, coh2))
                     , xlab = "Central Frequency (mHz)"
                     , ylab = "Offset Frequency (mHz)"
                     , legend.lab = expression(paste("Coherence (", sigma, " from mean)")))
  dev.off()

  png("~/school_lab/bin/testing/transfer3/plots/afterDecorrelate2.png"
      , height = 6, width = 8, res = 200, units = "in")
  par(mar = c(4,4,1,1))
  fields::image.plot(x = 1:dim(coh)[2], y = (-maxFreqOffsetIdx):maxFreqOffsetIdx, z=t(coh[,,1])
                     , zlim = c(-4,4)
                     , xlab = "Central Frequency Index", ylab = "Offset Frequency Index")
  dev.off()

  png("~/school_lab/bin/testing/transfer3/plots/afterDecorrelate3.png"
      , height = 6, width = 8, res = 200, units = "in")
  par(mar = c(4,4,1,1))
  fields::image.plot(x = 1:140, y = (-maxFreqOffsetIdx):maxFreqOffsetIdx, z=t(coh[,,1])
                     , zlim = c(-4,4)
                     , xlab = "Central Frequency Index", ylab = "Offset Frequency Index")
  dev.off()
}


hc.edm.decor <- function(){
  library(lubridate)
  library(fields)
  hcHomeDir <- "~/school_lab/contracts/health_canada2017/"
  mt24hr <- readRDS("~/school_lab/contracts/health_canada2017/assets/data/airPollution/mort24hrNoWeek-HighPass.rds")
  ap12hr <- readRDS("~/school_lab/contracts/health_canada2017/assets/data/airPollution/ap12hrNoWeek-LowPass-HighPass.rds")

  edmCd <- "4811"
  pols <- c("kO3.12h.lag0", "kNO2.12h.lag0", "kPM25A.12h.lag0")
  plotDir <- paste0(hcHomeDir, "assets/reports/13891694rttfkwhywqxp/images/")

  edm.tmp <- subset(ap12hr[[edmCd]], date <= ymd("2005-12-31") & date >= ymd("2005-01-01"))

  # standardize everything first:
  edm <- cbind(edm.tmp[, c("date", "amPm")], transfer2::df.std(edm.tmp[, pols]))
  mort <- subset(mt24hr[[edmCd]], date <= ymd("2005-12-31") & date >= ymd("2005-01-01"))

  x = edm[, pols]
  y = mort[, "Mort.CP.b.A0", drop = FALSE]


  nw = 8; k = 15; nFFTy = NULL; centre = "none"
  dty = 86400; dtx = 86400/2
  blockSizey = NULL; overlap = 0
  cohSigLev = 0.9; nOffAllowed = 3 ### careful here!
  forceZeroOffset = TRUE
  freqRange = (1e-6)*c(0, 5.8); maxFreqOffset = (1e-6)*1.8
  decorrelate = c("kPM25A.12h.lag0", "kNO2.12h.lag0", "kO3.12h.lag0")
  deadBand = c(0, 1/(58*86400)) # this is conservative
  forceZeroOffset = TRUE

  ## run the code in transfer3::tf up to, but not including, the decorrelation loop:
  # also need decorCentFreqIdx variable defined.
  # coh between 1 and 3 (O3 and PM)
  coh.b.o3pm <- transfer3::mscFromEigenHelper(yk1 = ykx.l$before[decorCentFreqIdx, , 1, 1]
                                            , yk2 = ykx.l$before[, , 1, 3]
                                            , k = k
                                            , nOffsets = nOffsets)
  coh.a.o3pm <- transfer3::mscFromEigenHelper(yk1 = ykx.l$after[decorCentFreqIdx, , 1, 1]
                                              , yk2 = ykx.l$after[, , 1, 3]
                                              , k = k
                                              , nOffsets = nOffsets)
  coh.a.o3pm[which(is.infinite(coh.a.o3pm))] <- min(coh.a.o3pm[!is.infinite(coh.a.o3pm)])

  png(paste0(plotDir, "plots/edmDecorO3Pm.png"), height = 5, width = 10, res = 201, units = "in")
  par(mar = c(4,4,1,2), mgp = c(2.5, 1, 0), mfrow = c(1, 2))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.b.o3pm), zlim = range(c(coh.b.o3pm, coh.a.o3pm))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.a.o3pm), zlim = range(c(coh.b.o3pm, coh.a.o3pm))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  dev.off()

  ###################################
  ## coh between 2 and 3 (NO2 and PM)
  ###################################
  coh.b.nopm <- transfer3::mscFromEigenHelper(yk1 = ykx.l$before[decorCentFreqIdx, , 1, 2]
                                              , yk2 = ykx.l$before[, , 1, 3]
                                              , k = k
                                              , nOffsets = nOffsets)
  coh.a.nopm <- transfer3::mscFromEigenHelper(yk1 = ykx.l$after[decorCentFreqIdx, , 1, 2]
                                              , yk2 = ykx.l$after[, , 1, 3]
                                              , k = k
                                              , nOffsets = nOffsets)
  coh.a.nopm[which(is.infinite(coh.a.nopm))] <- min(coh.a.nopm[!is.infinite(coh.a.nopm)])

  png(paste0(plotDir, "plots/edmDecorNoPm.png"), height = 5, width = 10, res = 201, units = "in")
  par(mar = c(4,4,1,2), mgp = c(2.5, 1, 0), mfrow = c(1, 2))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.b.nopm), zlim = range(c(coh.b.nopm, coh.a.nopm))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.a.nopm), zlim = range(c(coh.b.nopm, coh.a.nopm))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  dev.off()

  ###################################
  ## coh between 1 and 2 (O3 and NO2)
  ###################################
  coh.b.o3no <- transfer3::mscFromEigenHelper(yk1 = ykx.l$before[decorCentFreqIdx, , 1, 1]
                                              , yk2 = ykx.l$before[, , 1, 2]
                                              , k = k
                                              , nOffsets = nOffsets)
  coh.a.o3no <- transfer3::mscFromEigenHelper(yk1 = ykx.l$after[decorCentFreqIdx, , 1, 1]
                                              , yk2 = ykx.l$after[, , 1, 2]
                                              , k = k
                                              , nOffsets = nOffsets)
  coh.a.o3no[which(is.infinite(coh.a.o3no))] <- min(coh.a.o3no[!is.infinite(coh.a.o3no)])

  png(paste0(plotDir, "plots/edmDecorO3No.png"), height = 5, width = 10, res = 201, units = "in")
  par(mar = c(4,4,1,2), mgp = c(2.5, 1, 0), mfrow = c(1, 2))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.b.o3no), zlim = range(c(coh.b.o3no, coh.a.o3no))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  fields::image.plot(x = (1e6)*ykx.l$cFreq, y = (1e6)*ykx.l$oFreq
                     , z = t(coh.a.o3no), zlim = range(c(coh.b.o3no, coh.a.o3no))
                     , xlab = expression(paste("Centre Frequency (", mu, "Hz)"))
                     , ylab = expression(paste("Offset Frequency (", mu, "Hz)")))
  dev.off()
}

test.invertEigen <- function(){
  x = arima.sim(n = 200, model = list(ar = c(0.9, -0.2)))
  x.s <- spec.mtm(x, deltat = 2, dtUnits = "second", nw = 5, k = 9, returnInternals = TRUE)
  plot(x, type='l')

  x.hat1 <- invertEigenCoef(yk = x.s$mtm$eigenCoefs, v = x.s$mtm$dpss$v, dt = 2, nFFT = x.s$mtm$nFFT, N = 200)
  x.hat2 <- invertEigenCoef(yk = x.s$mtm$eigenCoefs * sqrt(x.s$mtm$eigenCoefWt), v = x.s$mtm$dpss$v, dt = 2, nFFT = x.s$mtm$nFFT, N = 200)
}
