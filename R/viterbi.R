
#' Viterbi algorithm
#'
#' For a given model, reconstructs the most probable states sequence,
#' using the Viterbi algorithm.
#'
#' @param m An object \code{moveHMM}
#' @param newdata An object \code{moveData} (optional)
#'
#' @return The sequence of most probable states.
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' # reconstruction of states sequence
#' states <- viterbi(m)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export

viterbi <- function(m, newdata = NULL)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")
    if(!is.null(newdata) && !is.moveData(newdata))
        stop("'newdata' must be a moveData object (as output by prepData) or NULL")

    data <- if (!is.null(newdata)) newdata else m$data
    nbStates <- ncol(m$mle$stepPar)
    beta <- m$mle$beta
    delta <- m$mle$delta
    stepDist <- m$conditions$stepDist
    angleDist <- m$conditions$angleDist
    stepPar <- m$mle$stepPar
    anglePar <- m$mle$anglePar
    zeroInflation <- m$conditions$zeroInflation
    knownStates <- m$knownStates

    if(nbStates==1)
        stop("No states to decode (nbStates=1)")

    # identify covariates
    covs <- model.matrix(m$conditions$formula, data)

    probs <- allProbs(data,nbStates,stepDist,angleDist,stepPar,anglePar,zeroInflation,knownStates)
    trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))

    nbAnimals <- length(unique(data$ID))
    # aInd = list of indices of first observation for each animal
    aInd <- c(1, which(data$ID[-1] != data$ID[-nrow(data)]) + 1)

    allStates <- NULL
    for(zoo in 1:nbAnimals) {
        nbObs <- length(which(data$ID==unique(data$ID)[zoo])) # nb of observations for animal zoo
        obsInd <- which(!is.na(data$step) & !is.na(data$angle))

        if(zoo!=nbAnimals) {
            p <- probs[aInd[zoo]:(aInd[zoo+1]-1),]
            tm <- trMat[,,aInd[zoo]:(aInd[zoo+1]-1)]
        }
        else {
            p <- probs[aInd[zoo]:nrow(probs),]
            tm <- trMat[,,aInd[zoo]:nrow(probs)]
        }

        xi <- matrix(NA,nbObs,nbStates)
        foo <- delta*p[1,]
        xi[1,] <- foo/sum(foo)
        for(i in 2:nbObs) {
            foo <- apply(xi[i-1,]*tm[,,i],2,max)*p[i,]
            xi[i,] <- foo/sum(foo)
        }

        stSeq <- rep(NA,nbObs)
        stSeq[nbObs] <- which.max(xi[nbObs,])
        for(i in (nbObs-1):1)
            stSeq[i] <- which.max(tm[,stSeq[i+1],i+1]*xi[i,])

        allStates <- c(allStates,stSeq)
    }

    return(allStates)
}
