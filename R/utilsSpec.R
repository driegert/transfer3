#' Calculate weighted eigencoefficients
#'
#' Calculates weighted eigencoefficients for a single vector / block
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
posFreq <- function(nFFT, dt){
  seq(0, 1/(2*dt), by = 1/(dt*nFFT))
}
