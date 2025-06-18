
#' Backward log-probabilities
#'
#' Used in \code{\link{stateProbs}}.
#'
#' @param m A \code{moveHMM} object.
#'
#' @return The matrix of backward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' lb <- logBeta(m)
#' }

logBeta <- function(m)
{
    data <- m$data
    nbStates <- ncol(m$mle$stepPar)
    nbAnimals <- length(unique(data$ID))
    nbObs <- nrow(data)
    lbeta <- matrix(NA,nbObs,nbStates)

    covs <- model.matrix(m$conditions$formula, m$data)

    probs <- allProbs(data,nbStates,m$conditions$stepDist,m$conditions$angleDist,m$mle$stepPar,
                      m$mle$anglePar,m$conditions$zeroInflation,m$knownStates)
    trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))

    aInd <- NULL
    for(i in 1:nbAnimals){
        aInd <- c(aInd,max(which(data$ID==unique(data$ID)[i])))
    }

    for(i in nbObs:1) {

        if(any(i==aInd)){
            foo <- rep(1,nbStates)
            lscale <- 0
        } else {
            gamma <- trMat[,,i+1]
            foo <- gamma%*%(probs[i+1,]*foo)
        }
        lbeta[i,] <- log(foo)+lscale
        sumfoo <- sum(foo)
        foo <- foo/sumfoo
        lscale <- lscale+log(sumfoo)
    }

    return(lbeta)
}
