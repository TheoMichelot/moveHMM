
#' Viterbi algorithm
#'
#' Reconstructs the most probable states sequence.
#'
#' @param data An object moveData.
#' @param nbStates Number of states of the HMM.
#' @param beta Matrix of regression coefficients for the transition probability matrix.
#' @param delta Stationary distribution.
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values. Defaults to "NULL"
#' if the turning angles distributions is not estimated.
#' @param stepPar Vector of state-dependent step length distribution parameters.
#' @param anglePar Vector of state-dependent turning angle distribution parameters.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#' @param zeroInflation TRUE if the step length distribution is inflated in zero.
#'
#' @return The sequence of most probable states.
#' @examples
#' data <- example$data
#' mod <- example$mod
#'
#' # reconstruction of states sequence
#' states <- viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar,
#'                   mod$mle$anglePar,angleMean)

viterbi <- function(data,nbStates,beta,delta,stepDist=c("gamma","weibull","exp"),
                    angleDist=c("NULL","vm","wrpcauchy"),stepPar,anglePar=NULL,angleMean=NULL,
                    zeroInflation=FALSE)
{
  # check arguments
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)

  if(nbStates<0) stop("nbStates should be at least 1.")
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$step)) stop("Missing field(s) in data.")

  vStepPar <- c(t(stepPar))
  vAnglePar <- c(t(anglePar))
  par <- c(stepPar,anglePar)
  p <- parDef(stepDist,angleDist,nbStates,is.null(angleMean),zeroInflation)
  bounds <- p$bounds
  parSize <- p$parSize
  if(sum(parSize)*nbStates!=length(par))
    stop("Wrong number of parameters.")
  stepBounds <- bounds[1:(parSize[1]*nbStates),]
  angleBounds <- bounds[(parSize[1]*nbStates+1):nrow(bounds),]
  if(length(which(vStepPar<stepBounds[,1] | vStepPar>stepBounds[,2]))>0 |
       length(which(vAnglePar<angleBounds[,1] | vAnglePar>angleBounds[,2]))>0)
    stop("Check the parameters bounds.")
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

  if(!is.null(angleMean))
    anglePar <- rbind(angleMean,anglePar)

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
