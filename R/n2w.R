
#' Scaling function: natural to working parameters.
#'
#' Scales each parameter from its natural interval to the set of real numbers, to allow for
#' unconstrained optimization. Used during the optimization of the log-likelihood.
#'
#' @param par Vector of state-dependent distributions parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in \code{par}. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param beta Matrix of regression coefficients for the transition probabilities.
#' @param delta Initial distribution. Default: \code{NULL} ; if the initial distribution is not estimated.
#' @param nbStates The number of states of the HMM.
#' @param estAngleMean \code{TRUE} if the angle mean is estimated, \code{FALSE} otherwise.
#'
#' @return A vector of unconstrained parameters.
#'
#' @examples
#' \dontrun{
#' nbStates <- 3
#' par <- c(0.001,0.999,0.5,0.001,1500.3,7.1)
#' bounds <- matrix(c(0,1, # bounds for first parameter
#'                    0,1, # bounds for second parameter
#'                    0,1, # ...
#'                    0,Inf,
#'                    0,Inf,
#'                    0,Inf),
#'                  byrow=TRUE,ncol=2)
#' beta <- matrix(rnorm(18),ncol=6,nrow=3)
#' delta <- c(0.6,0.3,0.1)
#'
#' # vector of working parameters
#' wpar <- n2w(par=par,bounds=bounds,beta=beta,delta=delta,nbStates=nbStates,
#'            estAngleMean=FALSE)
#' }
#'
#' @importFrom boot logit

n2w <- function(par,bounds,beta,delta=NULL,nbStates,estAngleMean)
{
  if(length(which(par<bounds[,1] | par>bounds[,2]))>0)
    stop("Check the parameters bounds.")

  nbPar <- length(par)/nbStates
  wpar <- NULL
  for(i in 1:nbPar) {
    index <- (i-1)*nbStates+1
    a <- bounds[index,1]
    b <- bounds[index,2]
    p <- par[index:(index+nbStates-1)]

    if(is.finite(a) & is.finite(b)) { # [a,b] -> R
      p <- logit((p-a)/(b-a))
    }
    else if(is.infinite(a) & is.finite(b)) { # ]-Inf,b] -> R
      p <- -log(-p+b)
    }
    else if(is.finite(a) & is.infinite(b)) { # [a,Inf[ -> R
      p <- log(p-a)
    }

    wpar <- c(wpar,p)
  }

  if(estAngleMean) {
    # identify angle distribution parameters
    foo <- length(wpar)-nbStates+1
    angleMean <- wpar[(foo-nbStates):(foo-1)]
    kappa <- wpar[foo:length(wpar)]

    # compute the working parameters for the angle distribution
    x <- kappa*cos(angleMean)
    y <- kappa*sin(angleMean)

    wpar[(foo-nbStates):(foo-1)] <- x
    wpar[foo:length(wpar)] <- y
  }

  wbeta <- as.vector(beta) # if beta is NULL, wbeta is NULL as well
  wdelta <- log(delta[-1]/delta[1]) # if delta is NULL, wdelta is NULL as well
  return(c(wpar,wbeta,wdelta))
}
