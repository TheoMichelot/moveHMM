
#' Generic CI method
#' @param m Fitted model
#' @param alpha Range of the confidence interval.
CI <- function(m,alpha) UseMethod("CI") # define generic method CI

#' Confidence intervals
#'
#' Computes the confidence intervals of the step length and turning angle parameters,
#' as well as for the transition probabilities regression parameters.
#'
#' @method CI moveHMM
#'
#' @param m A \code{moveHMM} object
#' @param alpha Range of the confidence intervals. Default : 0.95 (i.e. 95\% CIs).
#'
#' @return A list of the following objects :
#' \item{inf}{Inferior bound of the confidence interval for the parameters of the step lengths
#' distribution, for the parameters of the turning angle distribution, and for the transition
#' probabilities regression parameters}
#' \item{sup}{Superior bound of the confidence interval for the parameters of the step lengths
#' distribution, for the parameters of the turning angle distribution, and for the transition
#' probabilities regression parameters}
#'
#' @examples
#' m <- ex$m # moveHMM object, as returned by fitHMM
#'
#' CI(m)

CI.moveHMM <- function(m,alpha=0.95)
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

  est <- c(m$mod$estimate[1:i1],m$mod$estimate[i2:i3])

  var <- c(var[1:i1],var[i2:i3])

  # define appropriate quantile
  quantSup <- qnorm(1-(1-alpha)/2)

  # compute inf and sup for working parameters
  winf <- est-quantSup*sqrt(var)
  wsup <- est+quantSup*sqrt(var)

  # compute inf and sup on natural scale
  inf <- w2n(winf,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)
  sup <- w2n(wsup,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)

  # group CIs for step parameters, angle parameters, and t.p. coefficients
  inf <- list(stepPar=inf$stepPar,anglePar=angleCI(m,alpha)$inf,beta=inf$beta)
  sup <- list(stepPar=sup$stepPar,anglePar=angleCI(m,alpha)$sup,beta=sup$beta)

  return(list(inf=inf,sup=sup))
}
