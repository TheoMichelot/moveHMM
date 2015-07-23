
#' Fit an HMM to the data
#'
#' @param nbStates Number of states of the HMM.
#' @param data An object moveData.
#' @param stepPar0 Vector of initial state-dependent step length distribution parameters.
#' @param anglePar0 Vector of initial state-dependent turning angle distribution parameters.
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
#' stepPar <- c(15,10,50,20,0.2,0.3)
#' anglePar <- c(pi,0,0.7,2)
#' data <- simData(2,nbStates,"gamma","vm",stepPar,anglePar,nbCovs=2,zeroInflation=TRUE)
#'
#' mu0 <- c(20,80)
#' sigma0 <- c(20,40)
#' z0 <- c(0.5,0.5)
#' kappa0 <- c(0.1,0.1)
#' stepPar0 <- c(mu0,sigma0,z0)
#' anglePar0 <- kappa0
#'
#' nbCovs <- ncol(data[[1]]$covs)
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' delta0 <- c(1,1)/2
#'
#' mod <- fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,"gamma","vm",c(pi,0),TRUE)

fitHMM <- function(nbStates,data,stepPar0,anglePar0,beta0,delta0,stepDist=c("gamma","weibull","exp"),
                   angleDist=c("NULL","vm","wrpcauchy"),angleMean=NULL,
                   zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data[[1]]$step)) stop("Missing field(s) in data.")

  par0 <- c(stepPar0,anglePar0)

  p <- parDef(stepDist,angleDist,nbStates,is.null(angleMean),zeroInflation)
  bounds <- p$bounds
  parSize <- p$parSize
  if(sum(parSize)*nbStates!=length(par0))
    stop("Wrong number of initial parameters.")
  stepBounds <- bounds[1:(parSize[1]*nbStates),]
  angleBounds <- bounds[(parSize[1]*nbStates+1):nrow(bounds),]
  if(length(which(stepPar<stepBounds[,1] | stepPar>stepBounds[,2]))>0 |
       length(which(anglePar<angleBounds[,1] | anglePar>angleBounds[,2]))>0)
    stop("Check the parameters bounds.")

  if(ncol(data[[1]]$covs)>0) nbCovs <- ncol(data[[1]]$covs)
  else nbCovs <- 0

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)

  mle <- nlm(nLogLike,wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,
             zeroInflation,print.level=2,iterlim=1000)

  par <- w2n(mle$estimate,bounds,parSize,nbStates,nbCovs)
  return(par)
}
