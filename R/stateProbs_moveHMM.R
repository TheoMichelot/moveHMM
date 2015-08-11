
stateProbs <- function(m) UseMethod("stateProbs") # define generic method stateProbs

#' State probabilities
#'
#' @param m A moveHMM object.
#'
#' @return The matrix of state probabilities.
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#'
#' sp <- stateProbs(m)

stateProbs.moveHMM <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  la <- logAlpha(m)
  lb <- logBeta(m)
  c <- max(la[nbObs,])
  llk <- c + log(sum(exp(la[nbObs,]-c)))
  stateProbs <- matrix(NA,nbObs,nbStates)

  for(i in 1:nbObs)
    stateProbs[i,] <- exp(la[i,]+lb[i,]-llk)

  return(stateProbs)
}
