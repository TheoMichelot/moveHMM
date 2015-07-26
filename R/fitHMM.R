
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
#' # simulate data
#' nbAnimals <- 4
#' nbStates <- 2
#' nbCovs <- 3
#' mu<-c(20,70)
#' sigma<-c(8,20)
#' angleMean <- c(pi,0)
#' kappa <- c(0.8,1.3)
#' stepPar <- c(mu,sigma)
#' anglePar <- c(angleMean,kappa)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' zeroInflation <- FALSE
#'
#' data <- simData(nbAnimals,nbStates,stepDist,angleDist,stepPar,anglePar,nbCovs,zeroInflation)
#'
#' # estimation
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' angleMean <- c(pi,0)
#' stepPar0 <- c(mu0,sigma0)
#' anglePar0 <- kappa0
#'
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' delta0 <- rep(1,nbStates)/nbStates
#'
#' mod <- fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,"gamma","vm",angleMean,zeroInflation)

fitHMM <- function(nbStates,data,stepPar0,anglePar0,beta0,delta0,stepDist=c("gamma","weibull","exp"),
                   angleDist=c("NULL","vm","wrpcauchy"),angleMean=NULL,
                   zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field(s) in data.")

  par0 <- c(stepPar0,anglePar0)
  p <- parDef(stepDist,angleDist,nbStates,is.null(angleMean),zeroInflation)
  bounds <- p$bounds
  parSize <- p$parSize
  if(sum(parSize)*nbStates!=length(par0))
    stop("Wrong number of initial parameters.")
  stepBounds <- bounds[1:(parSize[1]*nbStates),]
  angleBounds <- bounds[(parSize[1]*nbStates+1):nrow(bounds),]
  if(length(which(stepPar0<stepBounds[,1] | stepPar0>stepBounds[,2]))>0 |
       length(which(anglePar0<angleBounds[,1] | anglePar0>angleBounds[,2]))>0)
    stop("Check the parameters bounds.")
  if(!is.null(angleBounds) & length(angleMean)!=nbStates)
    stop("The angleMean argument should be of length nbStates.")

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  nbCovs <- length(covsCol)

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)

  mle <- nlm(nLogLike,wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,
             zeroInflation,print.level=2,iterlim=1000)

  par <- w2n(mle$estimate,bounds,parSize,nbStates,nbCovs)
  return(par)
}
