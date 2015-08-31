
#' Generic deltaMethod method
#' @param m Fitted model
deltaMethod <- function(m) UseMethod("deltaMethod") # define generic method deltaMethod

#' Delta method
#'
#' Computes confidence intervals for the turning angle parameters, using the delta method.
#'
#' @method deltaMethod moveHMM
#'
#' @param m A \code{moveHMM} object.
#'
#' @return A list of :
#' \code{inf} inferior bound of 95% confidence interval;
#' \code{sup} superior bound of 95% confidence interval;
#' both for the parameters of the turning angle distribution.
#'
#' @examples
#' m <- ex$mod # moveHMM object, as returned by fitHMM
#'
#' deltaMethod(m)

deltaMethod.moveHMM <- function(m)
{
  if(m$angleDist=="none")
    stop("No angle parameter estimated.")

  if(length(m$mod)<=1)
    stop("The given model hasn't been fitted.")

  wpar <- m$mod$estimate
  nbStates <- ncol(m$mle$stepPar)

  # identify covariates
  covsCol <- which(names(m$data)!="ID" & names(m$data)!="x" & names(m$data)!="y" &
                     names(m$data)!="step" & names(m$data)!="angle")
  nbCovs <- length(covsCol)-1 # substract intercept column

  foo <- length(wpar)-nbStates+2
  wpar <- wpar[-(foo:length(wpar))]
  foo <- length(wpar)-(nbCovs+1)*nbStates*(nbStates-1)+1
  wpar <- wpar[-(foo:length(wpar))]
  foo <- length(wpar)-nbStates+1
  x <- wpar[(foo-nbStates):(foo-1)]
  y <- wpar[foo:length(wpar)]

  angleMean <- m$mle$anglePar[1,]
  kappa <- m$mle$anglePar[2,]

  Sigma <- ginv(m$mod$hessian)
  Sigma <- Sigma[(foo-nbStates):length(wpar),(foo-nbStates):length(wpar)]

  grad_kappa <- c(x/sqrt(x^2+y^2),y/sqrt(x^2+y^2))
  grad_angleMean <- c(-y/(x^2+y^2),1/(x+y^2/x))

  var_kappa <- grad_kappa%*%Sigma%*%grad_kappa
  if(m$conditions$estAngleMean)
    var_angleMean <- grad_angleMean%*%Sigma%*%grad_angleMean
  else
    var_angleMean <- 0

  inf <- rbind(angleMean-1.96*sqrt(var_angleMean),kappa-1.96*sqrt(var_kappa))
  sup <- rbind(angleMean+1.96*sqrt(var_angleMean),kappa+1.96*sqrt(var_kappa))

  return(list(inf=inf,sup=sup))
}
