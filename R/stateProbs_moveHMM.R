
#' State probabilities
#'
#' @param m A moveHMM object.
#'
#' @return The matrix of state probabilities.
#' @examples
#' m <- example$mod # moveHMM object (returned by fitHMM)
#'
#' sp <- stateProbs(m)

stateProbs <- function(m) UseMethod("stateProbs") # define generic method stateProbs

stateProbs.moveHMM <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  la <- lalpha(m)
  lb <- lbeta(m)
  c <- max(la[nbObs,])
  llk <- c + log(sum(exp(la[nbObs,]-c)))
  stateProbs <- matrix(NA,nbObs,nbStates)

  for(i in 1:nbObs) {
    cat(la[i,]," ; ",lb[i,]," ; ",llk,"\n")
    stateProbs[i,] <- exp(la[i,]+lb[i,]-llk)
  }

  return(stateProbs)
}
