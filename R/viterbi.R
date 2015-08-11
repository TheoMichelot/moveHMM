
#' Viterbi algorithm
#'
#' Reconstructs the most probable states sequence.
#'
#' @param data An object moveData.
#' @param nbStates Number of states of the HMM.
#' @param beta Matrix of regression coefficients for the transition probability matrix.
#' @param delta Stationary distribution.
#' @param stepDist Name of the distribution of the step lengths.
#' @param angleDist Name of the distribution of the turning angles.
#' Set to "none" if the angle distribution should not be estimated.
#' @param stepPar Vector of state-dependent step length distribution parameters.
#' @param anglePar Vector of state-dependent turning angle distribution parameters.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The sequence of most probable states.
#' @examples
#' data <- ex$data
#' m <- ex$mod
#' simPar <- ex$simPar
#'
#' # reconstruction of states sequence
#' states <- viterbi(data,simPar$nbStates,m$mle$beta,m$mle$delta,simPar$stepDist,
#'                   simPar$angleDist,m$mle$stepPar,m$mle$anglePar,simPar$angleMean)

viterbi <- function(data,nbStates,beta,delta,stepDist,angleDist,stepPar,
                    anglePar=NULL,angleMean=NULL,zeroInflation=FALSE)
{
  # check arguments
  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field(s) in data.")

  if(is.matrix(stepPar)) vStepPar <- c(t(stepPar))
  else vStepPar <- stepPar
  if(is.matrix(anglePar)) vAnglePar <- c(t(anglePar))
  else vAnglePar <- anglePar
  par <- c(stepPar,anglePar)
  p <- parDef(stepDist,angleDist,nbStates,TRUE,zeroInflation)
  bounds <- p$bounds
  parSize <- p$parSize

  if(sum(parSize)*nbStates!=length(par))
    stop("Wrong number of parameters.")

  stepBounds <- bounds[1:(parSize[1]*nbStates),]
  if(length(which(vStepPar<stepBounds[,1] | vStepPar>stepBounds[,2]))>0)
    stop("Check the step parameters bounds.")

  if(angleDist!="none") {
    angleBounds <- bounds[(parSize[1]*nbStates+1):nrow(bounds),]
    if(length(which(vAnglePar<angleBounds[,1] | vAnglePar>angleBounds[,2]))>0)
      stop("Check the angle parameters bounds.")
  }

  if(!is.null(angleMean) & length(angleMean)!=nbStates)
    stop("The angleMean argument should be of length nbStates.")

  # identify covariates
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  nbCovs <- length(covsCol)-1 # substract intercept column

  if(length(which(names(data)=="(Intercept)"))==0) { # no intercept column, if not called from fitHMM
    data <- cbind(data[,-covsCol],Intercept=rep(1,nrow(data)),data[,covsCol])
    covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                       names(data)!="step" & names(data)!="angle")
    nbCovs <- length(covsCol)-1 # substract intercept column
  }
  covs <- data[,covsCol]

  allProbs <- allProbs(data,nbStates,stepDist,angleDist,stepPar,anglePar,zeroInflation)
  trMat <- trMatrix_rcpp(nbStates,beta,as.matrix(covs))

  nbAnimals <- length(unique(data$ID))
  aInd <- NULL
  for(i in 1:nbAnimals)
    aInd <- c(aInd,which(data$ID==unique(data$ID)[i])[1])

  allStates <- NULL
  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(data$ID==unique(data$ID)[zoo])) # nb of observations for animal i
    obsInd <- which(!is.na(data$step) & !is.na(data$angle))

    if(zoo!=nbAnimals) {
      p <- allProbs[aInd[zoo]:(aInd[zoo+1]-1),]
      tm <- trMat[,,aInd[zoo]:(aInd[zoo+1]-1)]
    }
    else {
      p <- allProbs[aInd[zoo]:nrow(allProbs),]
      tm <- trMat[,,aInd[zoo]:nrow(allProbs)]
    }

    xi <- matrix(NA,nbObs,nbStates)
    foo <- delta*p[1,]
    xi[1,] <- foo/sum(foo)
    for(i in 2:nbObs) {
      foo <- apply(xi[i-1,]*tm[,,i],2,max)*p[i,]
      xi[i,] <- foo/sum(foo)
    }

    stSeq <- rep(NA,nbObs)
    stSeq[nbObs] <- which.max(xi[nbObs,])
    for(i in (nbObs-1):1)
      stSeq[i] <- which.max(tm[,stSeq[i+1],i]*xi[i,])

    allStates <- c(allStates,stSeq)
  }

  return(allStates)
}
