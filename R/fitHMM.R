
#' Fit an HMM to the data
#'
#' @param nbStates Number of states of the HMM.
#' @param data An object moveData.
#' @param stepPar0 Vector of initial state-dependent step length distribution parameters.
#' @param anglePar0 Vector of initial state-dependent turning angle distribution parameters.
#' @param beta0 Initial matrix of regression coefficients for the transition probability matrix.
#' @param delta0 Initial stationary distribution.
#' @param formula Regression formula for the covariates. Default : ~1 (no covariate).
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param angleMean Vector of state-dependent turning angles means. It defaults to NULL,
#' i.e. the means should be estimated.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The MLE of the parameters of the model.
#' @examples
#' data <- example$data
#' simPar <- example$simPar
#' par0 <- example$par0
#'
#' mod <- fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
#'               par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
#'               simPar$zeroInflation)

fitHMM <- function(nbStates,data,stepPar0,anglePar0,beta0=NULL,delta0=NULL,formula=~1,
                   stepDist=c("gamma","weibull","exp"),angleDist=c("NULL","vm","wrpcauchy"),
                   angleMean=NULL,zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field(s) in data.")

  par0 <- c(stepPar0,angleMean,anglePar0)
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
  if(!is.null(angleMean) & length(angleMean)!=nbStates)
    stop("The angleMean argument should be of length nbStates.")

  # build design matrix
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- model.matrix(formula,data)

  if(length(covsCol)>0) data <- cbind(data[-covsCol],covs)
  else data <- cbind(data,covs)
  nbCovs <- ncol(covs)-1 # substract intercept column

  # generate initial values for beta and delta
  if(is.null(beta0))
    beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                             nrow=nbCovs+1,byrow=TRUE)

  if(is.null(delta0)) delta0 <- rep(1,nbStates)/nbStates

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)

  # call to optimizer
  mod <- nlm(nLogLike,wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,
             zeroInflation,print.level=2,iterlim=1000)

  mle <- w2n(mod$estimate,bounds,parSize,nbStates,nbCovs)

  states <- viterbi(data,nbStates,mle$beta,mle$delta,stepDist,angleDist,mle$stepPar,mle$anglePar,
                 angleMean)

  mh <- list(data=data,mle=mle,states=states,stepDist=stepDist,angleDist=angleDist,
             angleMean=angleMean,mod=mod)
  return(moveHMM(mh))
}
