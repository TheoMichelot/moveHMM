
#' Stationary state probabilities
#'
#' Calculates the stationary probabilities of each state, for given
#' covariate values.
#'
#' @param m Fitted model (as output by \code{\link{fitHMM}}).
#' @param covMat Design matrix of covariates.
#'
#' @return Matrix of stationary state probabilities. Each row corresponds to
#' a row of covMat, and each column corresponds to a state.
#'
#' @export
stationary <- function(m, covMat)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    nbStates <- ncol(m$mle$stepPar)
    beta <- m$mle$beta

    if(nbStates==1)
        stop("No state probabilities (1-state model).")

    # all transition matrices
    allMat <- trMatrix_rcpp(nbStates=nbStates, beta=beta, covs=covMat)

    tryCatch({
        # for each transition matrix, derive corresponding stationary distribution
        probs <- apply(allMat, 3,
                       function(gamma)
                           solve(t(diag(nbStates)-gamma+1),rep(1,nbStates)))
        probs <- t(probs)
    },
    error = function(e) {
        stop(paste("The stationary probabilities cannot be calculated",
                   "for these covariate values (singular system)."))
    })

    return(probs)
}
