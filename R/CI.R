
#' Generic CI method
#' @param m Fitted model
#' @param alpha Range of the confidence interval.
#' @param nbSims Number of simulations.
#' @export
CI <- function(m,alpha,nbSims) UseMethod("CI")

#' Confidence intervals
#'
#' Computes the confidence intervals of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters.
#'
#' @method CI moveHMM
#'
#' @param m A \code{moveHMM} object
#' @param alpha Range of the confidence intervals. Default: 0.95 (i.e. 95\% CIs).
#' @param nbSims Number of simulations in the computation of the CIs for the angle parameters.
#' Default: 10^6.
#'
#' @return A list of the following objects:
#' \item{stepPar}{Confidence intervals for the parameters of the step lengths distribution}
#' \item{anglePar}{Confidence intervals for the parameters of the turning angles distribution}
#' \item{beta}{Confidence intervals for the regression coefficients of the transition probabilities.}
#'
#' @examples
#' m <- ex$m # moveHMM object, as returned by fitHMM
#'
#' CI(m)
#'
#' @export
#'
#' @importFrom MASS ginv

CI.moveHMM <- function(m,alpha=0.95,nbSims=10^6)
{
  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1.")

  nbStates <- ncol(m$mle$stepPar)

  # identify covariates
  covsCol <- which(names(m$data)!="ID" & names(m$data)!="x" & names(m$data)!="y" &
                     names(m$data)!="step" & names(m$data)!="angle")
  nbCovs <- length(covsCol)-1 # substract intercept column

  # inverse of Hessian
  Sigma <- ginv(m$mod$hessian)
  var <- diag(Sigma)

  p <- parDef(m$stepDist,m$angleDist,nbStates,m$conditions$estAngleMean,m$conditions$zeroInflation)

  # identify parameters of interest
  i1 <- p$parSize[1]*nbStates
  i2 <- sum(p$parSize)*nbStates+1
  i3 <- i2+nbStates*(nbStates-1)*(nbCovs+1)-1

  if(nbStates>1) {
    # select step parameters and "beta" parameters
    est <- c(m$mod$estimate[1:i1],m$mod$estimate[i2:i3])
    var <- c(var[1:i1],var[i2:i3])
  } else {
    # only select step parameters
    est <- m$mod$estimate[1:i1]
    var <- var[1:i1]
  }

  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)

  # compute lower and upper for working parameters
  wlower <- est-quantSup*sqrt(var)
  wupper <- est+quantSup*sqrt(var)

  # compute lower and upper on natural scale
  lower <- w2n(wlower,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)
  upper <- w2n(wupper,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)

  # CIs for angle parameters
  anglePar <- angleCI(m,alpha,nbSims)

  # group CIs for step parameters and t.p. coefficients
  stepPar <- list(lower=lower$stepPar,upper=upper$stepPar)
  beta <- list(lower=lower$beta,upper=upper$beta)

  if(!is.null(m$mle$beta))
    return(list(stepPar=stepPar,anglePar=anglePar,beta=beta))
  else
    return(list(stepPar=stepPar,anglePar=anglePar))
}
