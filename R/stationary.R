
#' Stationary state probabilities
#'
#' Calculates the stationary probabilities of each state, for given
#' covariate values.
#'
#' @param m Fitted model (as output by \code{\link{fitHMM}}).
#' @param covs Either a data frame or a design matrix of covariates.
#' @param beta Optional matrix of regression coefficients for the transition
#' probability model. By default, uses estimates in \code{m}.
#'
#' @return Matrix of stationary state probabilities. Each row corresponds to
#' a row of covs, and each column corresponds to a state.
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' # data frame of covariates
#' stationary(m, covs = data.frame(cov1 = 0, cov2 = 0))
#'
#' # design matrix (each column corresponds to row of m$mle$beta)
#' stationary(m, covs = matrix(c(1,0,cos(0)),1,3))
#'
#' @export
stationary <- function(m, covs, beta = m$mle$beta)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    nbStates <- ncol(m$mle$stepPar)

    if(nbStates==1)
        stop("No state probabilities (1-state model).")

    if(is.data.frame(covs)){
        covMat <- model.matrix(m$conditions$formula,covs)
    } else if(is.matrix(covs)){
        covMat <- covs
    } else stop("covs must either be a data frame or a matrix")

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
