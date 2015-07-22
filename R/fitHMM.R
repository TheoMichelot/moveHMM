
#' Fit an HMM to the data
#'
#' @param nbStates Number of states of the HMM.
#' @param data An object moveData.
#' @param par0 Vector of initial state-dependent distributions parameters.
#' @param beta0 Initial matrix of regression coefficients for the transition probability matrix.
#' @param delta0 Initial stationary distribution.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param angleMean Vector of state-dependent turning angles means. It defaults to NULL,
#' i.e. the means should be estimated.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The MLE of the parameters of the model.
#' @examples
#' nbStates <- 2
#' stepPar <- matrix(c(15,50,
#'                     10,20),
#'                   byrow=TRUE,ncol=2)
#' anglePar <- matrix(c(pi,0,
#'                      0.7,2),
#'                    byrow=TRUE,ncol=2)
#' data <- simData(2,nbStates,"gamma","vm",stepPar,anglePar,zeroInflation=0.3,nbCov=2)
#'
#' mu0 <- c(20,80)
#' sigma0 <- c(20,40)
#' kappa0 <- c(0.1,0.1)
#' z0 <- c(0.5,0.5)
#' par0 <- c(mu0,sigma0,z0,kappa0)
#'
#' nbCovs <- ncol(data[[1]]$covs)
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' delta0 <- c(1,1)/2
#'
#' mod <- fitHMM(nbStates,data,par0,beta0,delta0,"gamma","vm",c(pi,0),TRUE)

fitHMM <- function(nbStates,data,par0,beta0,delta0,stepDist=c("gamma","weibull","exp"),
                   angleDist=c("NULL","vm","wrpcauchy"),angleMean=NULL,
                   zeroInflation=FALSE)
{
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)

  parDef <- parDef(stepDist,angleDist,nbStates,is.null(angleMean),zeroInflation)
  parSize <- parDef$parSize
  if(sum(parSize)*nbStates!=length(par0))
    stop("Wrong number of initial parameters.")
  bounds <- parDef$bounds
  nbCovs <- ncol(data[[1]]$covs)

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)

  mle <- nlm(nLogLike,wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,
             zeroInflation,print.level=2,iterlim=1000)

  par <- w2n(mle$estimate,bounds,parSize,nbStates,nbCovs)
  return(par)
}
