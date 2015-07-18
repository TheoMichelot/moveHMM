
#' Scaling function : working to natural parameters
#'
#' Scales each parameter from the set of real numbers, back to its natural interval.
#' Used during the optimization of the log-likelihood.
#'
#' @param wpar Vector of state-dependent distributions  unconstrained parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in par. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param nbStates The number of states of the HMM.
#' @param nbCovs The number of covariates.
#'
#' @return A list containing a vector of natural (constrained) parameters, as well as delta
#' and beta.
#' @examples
#' nbStates <- 3
#' nbCovs <- 2
#' par <- c(0.001,0.999,0.5,0.001,1500.3,7.1)
#' bounds <- matrix(c(0,1,0,1,0,1,
#'                    0,Inf,0,Inf,0,Inf),
#'                  byrow=TRUE,ncol=2)
#' beta <- matrix(rnorm(18),ncol=6,nrow=3)
#' delta <- c(0.6,0.3,0.1)
#' wpar <- n2w(par,bounds,beta,delta,nbStates)
#' print(w2n(wpar,bounds,nbStates,nbCovs))
w2n <- function(wpar,bounds,nbStates,nbCovs)
{
  foo <- length(wpar)-nbStates+2
  delta <- wpar[foo:length(wpar)]
  delta <- exp(c(0,delta))
  delta <- delta/sum(delta)
  wpar <- wpar[-(foo:length(wpar))]

  foo <- length(wpar)-(nbCovs+1)*nbStates*(nbStates-1)+1
  beta <- wpar[foo:length(wpar)]
  beta <- matrix(beta,nrow=nbCovs+1)
  wpar <- wpar[-(foo:length(wpar))]

  nbPar <- length(wpar)/nbStates
  par <- NULL
  for(i in 1:nbPar) {
    index <- (i-1)*nbStates+1
    a <- bounds[index,1]
    b <- bounds[index,2]
    p <- wpar[index:(index+nbStates-1)]

    if(is.finite(a) & is.finite(b)) { # R -> [a,b]
      p <- (b-a)*inv.logit(p)+a
    }
    else if(is.infinite(a) & is.finite(b)) { # R -> ]-Inf,b]
      p <- exp(-p)-b
    }
    else if(is.finite(a) & is.infinite(b)) { # R -> [a,Inf[
      p <- exp(p)+a
    }

    par <- c(par,p)
  }

  return(list(par=par,beta=beta,delta=delta))
}
