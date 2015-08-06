
#' Decode the model
#'
#' @param m moveHMM object
decode.moveHMM <- function(m)
{
  data <- m$data
  mle <- m$mle
  stepDist <- m$stepDist
  angleDist <- m$angleDist
  nbStates <- ncol(mle$stepPar)
  angleMean <- m$angleMean
  states <- viterbi(data,nbStates,mle$beta,mle$delta,stepDist,angleDist,mle$stepPar,mle$anglePar,
                    angleMean)

  d <- list(states=states)
  return(d)
}
