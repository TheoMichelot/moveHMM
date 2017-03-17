
#' Forward log-probabilities
#'
#' Used in \code{\link{stateProbs}} and \code{\link{pseudoRes}}.
#'
#' @param m A \code{moveHMM} object.
#'
#' @return The matrix of forward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' la <- logAlpha(m)
#' }

logAlpha <- function(m)
{
    data <- m$data
    nbStates <- ncol(m$mle$stepPar)
    nbAnimals <- length(unique(data$ID))
    nbObs <- nrow(data)
    lalpha <- matrix(NA,nbObs,nbStates)

    beta <- m$mle$beta
    delta <- m$mle$delta

    covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                         names(data)!="step" & names(data)!="angle")
    covs <- data[,covsCol]

    probs <- allProbs(data,nbStates,m$conditions$stepDist,m$conditions$angleDist,m$mle$stepPar,
                      m$mle$anglePar,m$conditions$zeroInflation,m$knownStates)

    if(nbStates>1)
        trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))
    else
        trMat <- array(1,dim=c(1,1,nbObs))

    aInd <- NULL
    for(i in 1:nbAnimals)
        aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

    k <- 1
    for(i in 1:nbObs) {
        gamma <- trMat[,,i]
        if(any(i==aInd)){
            k <- max(nbAnimals,k+1)
            foo <- delta*probs[i,]
            lscale <- 0
        } else {
            foo <- (foo %*% gamma)*probs[i,]
        }
        lscale <- lscale+log(sum(foo))
        foo <- foo/sum(foo)
        lalpha[i,] <- log(foo)+lscale
    }

    return(lalpha)
}
