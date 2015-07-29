
#' Negative log-likelihood function
#'
#' @param wpar Vector of working parameters.
#' @param nbStates Number of states of the HMM.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in wpar. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values : c(number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution).
#' @param data An object moveData.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The negative log-likelihood of wpar given data.
#'
#' @examples
#' data <- example$data
#'
#' nbStates <- 2
#' nbCovs <- 2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' mu0 <- c(20,80)
#' sigma0 <- c(20,40)
#' kappa0 <- c(1,1)
#' beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
#'                 nrow=nbCovs+1,byrow=TRUE)
#' delta0 <- c(1,1)/2
#' par0 <- c(mu0,sigma0,kappa0)
#' wpar <- n2w(par0,bounds,beta0,delta0,nbStates)
#' angleMean <- c(pi,0)
#'
#' bounds <- parDef(stepDist,angleDist,nbStates,FALSE,FALSE)$bounds
#' parSize <- parDef(stepDist,angleDist,nbStates,FALSE,FALSE)$parSize
#'
#' l <- nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,FALSE)

nLogLike <- function(wpar,nbStates,bounds,parSize,data,stepDist=c("gamma","weibull","exp"),
                     angleDist=c("NULL","vm","wrpcauchy"),angleMean=NULL,zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)
  if(nbStates<1) stop("nbStates must be at least 1.")
  if(length(wpar)!=sum(parSize)*nbStates+nbStates*(nbStates-1)*(nbCovs+1)+nbStates-1)
    stop("Wrong number of parameters in wpar.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field(s) in data.")

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  nbCovs <- length(covsCol)-1 # substract intercept column

  if(length(which(names(data)=="(Intercept)"))==0) { # no intercept column, if not called from fitHMM
    data <- cbind(data[,-covsCol],Intercept=rep(1,nrow(data)),data[,covsCol])
    covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                       names(data)!="step" & names(data)!="angle")
    nbCovs <- length(covsCol)-1 # substract intercept column
  }

  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs)

  if(!is.null(angleMean)) # if the turning angles' mean is not estimated
    par$anglePar <- rbind(angleMean,par$anglePar)

  nbObs <- length(data$step)
  covs <- data[,covsCol]
  trMat <- trMatrix_rcpp(nbStates,par$beta,as.matrix(covs))
  allProbs <- allProbs(data,nbStates,stepDist,angleDist,par$stepPar,par$anglePar,zeroInflation)

  nbAnimals <- length(unique(data$ID))
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

  llk <- 0
  for(i in 1:nbAnimals)
  {
    if(i!=nbAnimals) {
      p <- allProbs[aInd[i]:(aInd[i+1]-1),]
      tm <- trMat[,,aInd[i]:(aInd[i+1]-1)]
    }
    else {
      p <- allProbs[aInd[i]:nrow(allProbs),]
      tm <- trMat[,,aInd[i]:nrow(allProbs)]
    }
    lscale <- nLogLike_rcpp(tm,par$delta,p) # call to C++ function
    llk <- llk+lscale
  }
  return(-llk)
}
