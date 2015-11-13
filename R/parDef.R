
#' Parameters definition
#'
#' @param stepDist Name of the distribution of the step lengths.
#' @param angleDist Name of the distribution of the turning angles.
#' Set to \code{"none"} if the angle distribution should not be estimated.
#' @param nbStates Number of states of the HMM.
#' @param estAngleMean \code{TRUE} if the mean of the turning angles distribution is estimated,
#' \code{FALSE} otherwise.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#'
#' @return A list of:
#' \item{parSize}{Vector of two values: number of parameters of the step length distribution,
#' number of parameters of the turning angle distribution}
#' \item{bounds}{Matrix with 2 columns and \code{sum(parSize)} rows - each row contains the lower and upper
#' bound for the correponding parameter)}
#' \item{parNames}{Names of parameters of step distribution (the names of the parameters of the
#' angle distribution are always the same).}

parDef <- function(stepDist,angleDist,nbStates,estAngleMean,zeroInflation)
{
  parSize <- c(NA,NA)

  switch(stepDist,
         "gamma"={
           parSize[1] <- 2
           stepBounds <- matrix(rep(c(0,Inf),2*nbStates),ncol=2,byrow=TRUE)
           parNames <- c("mean","sd")
         },
         "weibull"={
           parSize[1] <- 2
           stepBounds <- matrix(rep(c(0,Inf),2*nbStates),ncol=2,byrow=TRUE)
           parNames <- c("shape","scale")
         },
         "lnorm"={
           parSize[1] <- 2
           stepBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(0,Inf),nbStates)),
                                ncol=2,byrow=TRUE)
           parNames <- c("location","scale")
         },
         "exp"={
           parSize[1] <- 1
           stepBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
           parNames <- c("rate")
         })

  # include zero-mass
  if(zeroInflation) {
    parSize[1] <- parSize[1]+1
    stepBounds <- rbind(stepBounds,matrix(rep(c(0,1),nbStates),ncol=2,byrow=TRUE))
    parNames <- c(parNames,"zero-mass")
  }

  switch(angleDist,
         "none"={
           parSize[2] <- 0
           angleBounds <- NULL
         },
         "vm"={
           if(estAngleMean) { # if the angle mean is estimated
             parSize[2] <- 2
             # bounds are chosen such that the parameters are not scaled
             # (already in the right intervals for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,Inf),nbStates)),
                                   ncol=2,byrow=TRUE)
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0,Inf),nbStates),ncol=2,byrow=TRUE)
           }
         },
         "wrpcauchy"={
           if(estAngleMean) {
             parSize[2] <- 2
             # bounds are chosen such that the mean is not scaled, but the concentration is
             # scaled from ]0,1[ to ]0,Inf[ (for computing x and y)
             angleBounds <- matrix(c(rep(c(-Inf,Inf),nbStates),rep(c(-Inf,1),nbStates)),
                                   ncol=2,byrow=TRUE)
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(rep(c(0,1),nbStates),ncol=2,byrow=TRUE)
           }
         })

  bounds <- rbind(stepBounds,angleBounds)
  return(list(parSize=parSize,bounds=bounds,parNames=parNames))
}
