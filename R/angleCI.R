
#' Confidence intervals for angle parameters
#'
#' Simulation-based computation of confidence intervals for the parameters of the angle distribution.
#' Used in \code{\link{CI.moveHMM}}.
#'
#' @param m A \code{moveHMM} object
#' @param alpha Range of the confidence intervals. Default : 0.95 (i.e. 95\% CIs).
#'
#' @return A list of the following objects :
#' \item{inf}{Inferior bound of the confidence interval for the parameters of the angle distribution}
#' \item{sup}{Superior bound of the confidence interval for the parameters of the angle distribution}

angleCI <- function(m,alpha)
{
  if(alpha<0 | alpha>1)
    stop("alpha needs to be between 0 and 1")

  nbStates <- ncol(m$mle$anglePar)
  inf <- matrix(NA,ncol=nbStates,nrow=2)
  sup <- matrix(NA,ncol=nbStates,nrow=2)

  pdef <- parDef(m$stepDist,m$angleDist,nbStates,estAngleMean=TRUE,m$conditions$zeroInflation)
  parSize <- pdef$parSize
  bounds <- pdef$bounds

  for(state in 1:nbStates) {
    nbSims <- 10^6 # the bigger the better

    # working MLE
    wpar <- m$mod$estimate
    # only keep the angle parameters
    wpar <- wpar[(parSize[1]*nbStates+1):(sum(parSize)*nbStates)]
    x <- wpar[1:(length(wpar)/2)]
    y <- wpar[(length(wpar)/2+1):length(wpar)]

    # compute cov matrix
    Sigma <- ginv(m$mod$hessian)
    # only keep row/columns for angle parameters
    Sigma <- Sigma[(parSize[1]*nbStates+1):(sum(parSize)*nbStates),
                   (parSize[1]*nbStates+1):(sum(parSize)*nbStates)]
    # only keep row/columns for current state
    ind <- c(state,state+nbStates)
    Sigma <- Sigma[ind,ind]

    # simulated working parameters
    wSims <- mvrnorm(nbSims, mu=c(x[state],y[state]), Sigma=Sigma)

    # simulated natural parameters
    nSims <- cbind(Arg(wSims[,1]+1i*wSims[,2]),
                   sqrt(wSims[,1]^2+wSims[,2]^2))

    # scale concentration if necessary
    # (assumes that concentration has bounds like ]-Inf,b])
    if(is.finite(bounds[sum(parSize)*nbStates,2]))
      nSims[,2] <- -(exp(-nSims[,2])-bounds[sum(parSize)*nbStates,2])

    # if the mean is around -pi/pi, add 2*pi to negative values
    if(m$mle$anglePar[1,state]<(-pi/2) | m$mle$anglePar[1,state]>(pi/2))
      nSims[which(nSims[,1]<0),1] <- nSims[which(nSims[,1]<0),1]+2*pi

    # define appropriate quantile
    quantInf <- (1-alpha)/2
    quantSup <- 1-quantInf

    # compute CIs
    inf[1,state] <- quantile(nSims[,1],quantInf)
    inf[2,state] <- quantile(nSims[,2],quantInf)
    sup[1,state] <- quantile(nSims[,1],quantSup)
    sup[2,state] <- quantile(nSims[,2],quantSup)
  }

  return(list(inf=inf,sup=sup))
}
