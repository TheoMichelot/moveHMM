
#' Backward log-probabilities
#'
#' Used in stateProbs.
#'
#' @param m A moveHMM object.
#'
#' @return The matrix of backward log-probabilities.
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#'
#' lb <- logBeta(m)

logBeta <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  lbeta <- matrix(NA,nbObs,nbStates)

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- data[,covsCol]

  allProbs <- allProbs(data,nbStates,m$stepDist,m$angleDist,m$mle$stepPar,m$mle$anglePar,m$zeroInflation)
  trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))

  lscale <- log(nbStates)
  foo <- rep(1,nbStates)/nbStates
  lbeta[nbObs,] <- rep(0,nbStates)

  for(i in (nbObs-1):1) {
    gamma <- trMat[,,(i+1)]
    foo <- gamma%*%(allProbs[i+1,]*foo)
    lbeta[i,] <- log(foo)+lscale
    sumfoo <- sum(foo)
    foo <- foo/sumfoo
    lscale <- lscale+log(sumfoo)
  }

  return(lbeta)
}
