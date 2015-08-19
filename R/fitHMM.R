
#' Fit an HMM to the data
#'
#' @param nbStates Number of states of the HMM.
#' @param data An object moveData.
#' @param stepPar0 Vector of initial state-dependent step length distribution parameters.
#' The parameters should be in the order expected by the pdf of stepDist, and the (optional) zero-mass
#' parameter should be the last. For example, for a 2-state model using the Gamma (gamma) distribution and
#' including zero-inflation, the vector of initial parameters would be something like :
#' c(mu1,mu2,sigma1,sigma2,zeromass1,zeromass2).
#' @param anglePar0 Vector of initial state-dependent turning angle distribution parameters.
#' The parameters should be in the order expected by the pdf of angleDist. For example, for a 2-state
#' model using the Von Mises (vm) distribution, the vector of initial parameters would be something like :
#' c(mu1,mu2,kappa1,kappa2).
#' @param beta0 Initial matrix of regression coefficients for the transition probabilities. Default : NULL.
#' If not specified, beta0 is initialized such that the diagonal elements of the transition probability
#' matrix are dominant.
#' @param delta0 Initial value for the initial distribution of the HMM. Default : rep(1/nbStates,nbStates).
#' @param formula Regression formula for the covariates. Default : ~1 (no covariate effect).
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are : gamma, weibull, lnorm, exp. Default : gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are : vm, wrpcauchy. Set to "none" if the angle distribution should
#' not be estimated. Default : vm.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Default : NULL (the angle mean is estimated).
#' @param zeroInflation TRUE if the step length distribution is inflated in zero. Default : FALSE. If TRUE,
#' initial values for the zero-mass parameters should be included in stepPar0.
#' @param stationary FALSE if there are covariates. If TRUE, the initial distribution is considered
#' equal to the stationary distribution. Default : FALSE.
#' @param verbose Determines the print level of the optimizer. The default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed.
#'
#' @return A moveHMM object, including the MLE of the model parameters.
#' @examples
#' ### 1. simulate data
#' # define all the arguments of simData
#' nbAnimals <- 2
#' nbStates <- 2
#' nbCovs <- 2
#' mu<-c(15,50)
#' sigma<-c(10,20)
#' angleMean <- c(pi,0)
#' kappa <- c(0.7,1.5)
#' stepPar <- c(mu,sigma)
#' anglePar <- c(angleMean,kappa)
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' zeroInflation <- FALSE
#' obsPerAnimal <- c(50,100)
#'
#' data <- simData(nbAnimals,nbStates,stepDist,angleDist,stepPar,anglePar,NULL,nbCovs,zeroInflation,
#'                 obsPerAnimal)
#'
#' ### 2. fit the model to the simulated data
#' # define initial values for the parameters
#' mu0 <- c(20,70)
#' sigma0 <- c(10,30)
#' kappa0 <- c(1,1)
#' stepPar0 <- c(mu0,sigma0) # no zero-inflation, so no zero-mass included
#' anglePar0 <- kappa0 # the angle mean is not estimated, so only the concentration parameter is needed
#' formula <- ~cov1+cos(cov2)
#'
#' mod <- fitHMM(nbStates,data,stepPar0,anglePar0,beta0=NULL,delta0=NULL,formula,
#'               stepDist="gamma",angleDist="vm",angleMean,zeroInflation,verbose=2)

fitHMM <- function(nbStates,data,stepPar0,anglePar0,beta0=NULL,delta0=NULL,formula=~1,
                   stepDist=c("gamma","weibull","lnorm","exp"),angleDist=c("vm","wrpcauchy","none"),
                   angleMean=NULL,zeroInflation=FALSE,stationary=FALSE,verbose=0)
{
  # build design matrix
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- model.matrix(formula,data)

  if(length(covsCol)>0) data <- cbind(data[-covsCol],covs)
  else data <- cbind(data,covs)
  nbCovs <- ncol(covs)-1 # substract intercept column

  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field in data : step.")

  par0 <- c(stepPar0,anglePar0)
  p <- parDef(stepDist,angleDist,nbStates,is.null(angleMean),zeroInflation)
  bounds <- p$bounds
  parSize <- p$parSize
  if(sum(parSize)*nbStates!=length(par0)) {
    error <- "Wrong number of initial parameters."
    if(parSize[1]*nbStates!=length(stepPar0))
      error <- paste(error,": there should be",parSize[1]*nbStates,"initial step parameters.")
    if(angleDist!="none" & parSize[2]*nbStates!=length(stepPar0))
      error <- paste(error,": there should be",parSize[2]*nbStates,"initial angle parameters.")
    stop(error)
  }

  if(!is.null(beta0)) {
    if(ncol(beta0)!=nbStates*(nbStates-1) | nrow(beta0)!=nbCovs+1) {
      error <- paste("beta0 has wrong dimensions : it should have",nbCovs+1,"rows and",
                     nbStates*(nbStates-1),"columns.")
      stop(error)
    }
  }

  if(!is.null(delta0))
    if(length(delta0)!=nbStates)
      stop(paste("delta0 has the wrong length : it should have",nbStates,"elements."))

  stepBounds <- bounds[1:(parSize[1]*nbStates),]
  if(length(which(stepPar0<stepBounds[,1] | stepPar0>stepBounds[,2]))>0)
    stop("Check the step parameters bounds.")

  if(angleDist!="none") {
    angleBounds <- bounds[(parSize[1]*nbStates+1):nrow(bounds),]
    if(length(which(anglePar0<angleBounds[,1] | anglePar0>angleBounds[,2]))>0)
      stop("Check the angle parameters bounds.")
  }
  if(!is.null(angleMean) & length(angleMean)!=nbStates)
    stop("The angleMean argument should be of length nbStates.")

  # check that observations are within expected bounds
  if(length(which(data$step<0))>0)
    stop("The step lengths should be positive.")
  if(length(which(data$angle < -pi | data$angle > pi))>0)
    stop("The turning angles should be between -pi and pi.")

  # check that zeroInflation is consistent with the observations
  if(!zeroInflation & length(which(data$step==0))>0)
    stop("Zero-inflation should be included if step length can be zero.")
  if(zeroInflation & length(which(data$step==0))==0)
    stop("Zero-inflation should not be included if step length is never zero.")

  # check that stationary==FALSE if there are covariates
  if(nbCovs>0 & stationary==TRUE)
    stop("stationary can't be set to TRUE if there are covariates.")

  # generate initial values for beta and delta
  if(is.null(beta0))
    beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                             nrow=nbCovs+1,byrow=TRUE)

  if(is.null(delta0)) delta0 <- rep(1,nbStates)/nbStates
  if(stationary) delta0 <- NULL

  estAngleMean <- (is.null(angleMean) & angleDist!="none")

  # build the vector of working parameters
  wpar <- n2w(par0,bounds,beta0,delta0,nbStates,estAngleMean)

  # this function is used to muffle the warning "NA/Inf replaced by maximum positive value" in nlm
  h <- function(w) {
    if(any(grepl("NA/Inf replaced by maximum positive value",w)))
      invokeRestart("muffleWarning")
  }

  # call to optimizer nlm
  withCallingHandlers(mod <- nlm(nLogLike,wpar,nbStates,bounds,parSize,data,stepDist,
                                 angleDist,angleMean,zeroInflation,stationary,
                                 print.level=verbose,
                                 iterlim=1000,
                                 hessian=TRUE),
                      warning=h) # filter warnings using function h

  # convert the parameters back to their natural scale
  mle <- w2n(mod$estimate,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary)

  if(stationary) {
    gamma <- trMatrix_rcpp(nbStates,mle$beta,covs)[,,1]
    mle$delta <- solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))
  }

  if(!is.null(angleMean) & angleDist!="none") {
    mle$anglePar <- rbind(angleMean,mle$anglePar)
    rownames(mle$anglePar) <- NULL # remove rbind row name
  }

  # decode the sequence of states
  states <- viterbi(data,nbStates,mle$beta,mle$delta,stepDist,angleDist,mle$stepPar,mle$anglePar,
                    angleMean,zeroInflation)

  mh <- list(data=data,mle=mle,stepDist=stepDist,angleDist=angleDist,
             mod=mod,states=states,zeroInflation=zeroInflation,
             estAngleMean=estAngleMean,stationary=stationary)
  return(moveHMM(mh))
}
