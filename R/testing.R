# testing the Fortran code to make sure it makes sense!

#' @useDynLib transfer3
test.cohMsc <- function(){
  set.seed(19840331)
  yk1 <- matrix(rnorm(12), ncol = 3, nrow = 4)
  yk2 <- matrix(rnorm(24), ncol = 3, nrow = 8)

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

}
