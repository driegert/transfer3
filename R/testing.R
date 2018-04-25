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
