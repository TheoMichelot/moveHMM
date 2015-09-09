
#' Generic stateProbs method
#' @param m Fitted model
stateProbs <- function(m) UseMethod("stateProbs") # define generic method stateProbs

#' State probabilities
#' @method stateProbs moveHMM
#'
#' @param m A \code{moveHMM} object.
#'
#' @return The matrix of state probabilities.
#'
#' @examples
#' m <- ex$m # moveHMM object (returned by fitHMM)
#'
#' sp <- stateProbs(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).

stateProbs.moveHMM <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  la <- logAlpha(m) # forward log-probabilities
  lb <- logBeta(m) # backward log-probabilities
  c <- max(la[nbObs,]) # cancels out below ; prevents numerical errors
  llk <- c + log(sum(exp(la[nbObs,]-c)))
  stateProbs <- matrix(NA,nbObs,nbStates)

  for(i in 1:nbObs)
    stateProbs[i,] <- exp(la[i,]+lb[i,]-llk)

  return(stateProbs)
}
