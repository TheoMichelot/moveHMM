
#' Fit an HMM to the data
#'
#' @param nbStates Number of states of the HMM.
#' @param data An object moveData.
#' @param stepPar0 Vector of initial state-dependent step length distribution parameters.
#' @param anglePar0 Vector of initial state-dependent turning angle distribution parameters.
#' @param beta0 Initial matrix of regression coefficients for the transition probability matrix.
#' @param delta0 Initial stationary distribution.
#' @param formula Regression formula for the covariates. Default : ~1 (no covariate).
#' @param stepDist Name of the distribution of the step lengths.
#' Supported distributions are : gamma, weibull, lnorm, exp.
#' @param angleDist Name of the distribution of the turning angles.
#' Supported distributions are : vm, wrpcauchy. Set to "none" if the angle distribution should
#' not be estimated.
#' @param angleMean Vector of state-dependent turning angles means. It defaults to NULL,
#' i.e. the means should be estimated.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#' @param stationary FALSE if there are covariates. If TRUE, the initial distribution is considered
#' equal to the stationary distribution.
#' @param verbose Determines the print level of the optimizer. The default value of 0 means that no
#' printing occurs, a value of 1 means that the first and last iterations of the optimization are
#' detailed, and a value of 2 means that each iteration of the optimization is detailed.
#'
#' @return The MLE of the parameters of the model.
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
#' stepPar0 <- c(mu0,sigma0)
#' anglePar0 <- kappa0
#' formula <- ~cov1+cos(cov2)
#'
#' mod <- fitHMM(nbStates,data,stepPar0,anglePar0,NULL,NULL,formula,
#'               "gamma","vm",angleMean,zeroInflation,verbose=2)

fitHMM <- function(nbStates,data,stepPar0,anglePar0,beta0=NULL,delta0=NULL,formula=~1,
                   stepDist=c("gamma","weibull","lnorm","exp"),angleDist=c("vm","wrpcauchy","none"),
                   angleMean=NULL,zeroInflation=FALSE,stationary=FALSE,verbose=0)
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

  # build design matrix
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- model.matrix(formula,data)

  if(length(covsCol)>0) data <- cbind(data[-covsCol],covs)
  else data <- cbind(data,covs)
  nbCovs <- ncol(covs)-1 # substract intercept column

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

  mle <- w2n(mod$estimate,bounds,parSize,nbStates,nbCovs,estAngleMean,stationary)

  if(stationary) {
    gamma <- trMatrix_rcpp(nbStates,mle$beta,covs)[,,1]
    mle$delta <- solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))
  }

  if(!is.null(angleMean) & angleDist!="none")
    mle$anglePar <- rbind(angleMean,mle$anglePar)

  states <- viterbi(data,nbStates,mle$beta,mle$delta,stepDist,angleDist,mle$stepPar,mle$anglePar,
                    angleMean,zeroInflation)

  mh <- list(data=data,mle=mle,stepDist=stepDist,angleDist=angleDist,
             mod=mod,states=states,zeroInflation=zeroInflation)
  return(moveHMM(mh))
}
