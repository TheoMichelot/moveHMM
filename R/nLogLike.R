
#' Negative log-likelihood function
#'
#' @param nbStates Number of states of the HMM.
#' @param wpar Vector of working parameters.
#' @param bounds Matrix with 2 columns and as many rows as there are elements in wpar. Each row
#' contains the lower and upper bound for the correponding parameter.
#' @param parSize Vector of two values : c(number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution).
#' @param data An object moveData.
#' @param angleMean Vector of means of turning angles if not estimated (one for each state).
#' Defaults to NULL.
#'
#' @return The negative log-likelihood of wpar given data.
nLogLike <- function(nbStates,wpar,bounds,parSize,data,angleMean=NULL)
{
  llk <- 0
  nbAnimals <- length(data)
  if(!is.null(data[[1]]$covs)) nbCovs <- ncol(data[[1]]$covs)
  else nbCovs <- 0

  par <- w2n(wpar,bounds,parSize,nbStates,nbCovs)
  if(!is.null(angleMean)) # if the turning angles' mean is not estimated
    par$anglePar <- rbind(angleMean,par$anglePar)

  for(zoo in 1:nbAnimals) {
    nbObs <- length(data[[zoo]]$step)
    covs <- data[[zoo]]$covs

    trMat <- trMatrix(nbStates,nbObs,par$beta,covs)
    allProbs <- allProbs(data[[zoo]],nbStates,stepFun,angleFun,par$stepPar,par$anglePar)

    lscale <- 0
    alpha <- par$delta*allProbs[1,]
    for(k in 2:nbObs) {
      gamma <- trMat[,,k]
      alpha <- alpha%*%gamma*allProbs[k,]

      # scaling
      sumalpha <- sum(alpha)
      lscale <- lscale+log(sumalpha)
      alpha <- alpha/sumalpha
    }
    llk <- llk + lscale
  }

  return(-llk)
}
