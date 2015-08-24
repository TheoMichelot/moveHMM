
#' Generic confIntervals method
#' @param m Fitted model
confIntervals <- function(m) UseMethod("confIntervals") # define generic method confIntervals

#' Confidence intervals
#'
#' Computes the confidence intervals of the step length parameters, as well as for the transition
#' probabilities regression parameters.
#'
#' @method confIntervals moveHMM
#'
#' @param m A \code{moveHMM} object
#'
#' @return A list of the following objects :
#' \item{inf}{Inferior bound of 95% confidence interval for the parameters of the step lengths
#' distribution, and for the transition probabilities regression parameters}
#' \item{sup}{Superior bound of 95% confidence interval for the parameters of the step lengths
#' distribution, and for the transition probabilities regression parameters}
#'
#' @examples
#' m <- ex$mod # moveHMM object, as returned by fitHMM
#'
#' confIntervals(m)

confIntervals.moveHMM <- function(m)
{
  if(is.na(m$mod))
    stop("The given model hasn't been fitted.")

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

  # compute inf and sup for working parameters
  winf <- est-1.96*sqrt(var)
  wsup <- est+1.96*sqrt(var)

  # compute inf and sup on natural scale
  inf <- w2n(winf,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)
  sup <- w2n(wsup,p$bounds[1:i1,],c(p$parSize[1],0),nbStates,nbCovs,FALSE,TRUE)

  inf <- list(stepPar=inf$stepPar,beta=inf$beta)
  sup <- list(stepPar=sup$stepPar,beta=sup$beta)

  return(list(inf=inf,sup=sup))
}
