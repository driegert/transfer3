# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#' @export
hello <- function() {
  print("Hello, world!")
}

# really... any of these methods are acceptable I suppose.  12 seconds is a lot of seconds though..
#' @export
testSpecMtmTime <- function(){
  n <- 11520
  nw <- 8
  k <- 15
  dt <- 16*60
  dat <- data.frame(matrix(rnorm(n*3*40), nrow = 11520))

  res <- list()
  sleps <- dpss(n = n, k = k, nw = nw, returnEigenvalues = FALSE)$v

  timeStart <- proc.time()
  # 1)
  for (i in 1:ncol(dat)){
    res[[i]] <- spec.mtm(dat[, i], centre = "none", plot = FALSE, returnInternals = TRUE
                         , nw = nw, k = k, deltat = dt, dtUnits = "second"
                         , dpssIN = sleps)$mtm$eigenCoefs
  }

  # 2)
  # res <- lapply(dat, spec.mtm, centre = "none", plot = FALSE, returnInternals = TRUE
  #                                      , nw = nw, k = k, deltat = dt, dtUnits = "second"
  #                                      , dpssIN = sleps)

  # 3)
  # cl <- makeCluster(rep("localhost", 4), type = "SOCK")
  # clusterSetupRNG(cl)
  # res <- clusterApply(cl, dat, spec.mtm, centre = "none", plot = FALSE, returnInternals = TRUE
  #                                                          , nw = nw, k = k, deltat = dt, dtUnits = "second"
  #                                                          , dpssIN = sleps)
  # stopCluster(cl)

  endTime <- proc.time()

  print(endTime - timeStart)
}
