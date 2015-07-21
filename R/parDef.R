
#' Parameters definition
#'
#' @param stepDist Name of the distribution of the step length values.
#' @param angleDist Name of the distribution of the turning angle values.
#' @param nbStates Number of states of the HMM.
#' @param estAngleMean TRUE if the mean of the turning angles distribution is estimated,
#' FALSE otherwise.
#'
#' @return A list of two elements : parSize (= c(number of parameters of the step length
#' distribution,number of parameters of the turning angle distribution)) and bounds
#' (matrix with 2 columns and sum(parSize) rows. Each row contains the lower and upper
#' bound for the correponding parameter).

parDef <- function(stepDist=c("gamma","weibull","exp"),angleDist=c("vm","wrpcauchy"),nbStates,
                   estAngleMean)
{
  stepDist <- match.arg(stepDist)
  angleDist <- match.arg(angleDist)

  parSize <- c(NA,NA)
  switch(stepDist,
         "gamma"={
           parSize[1] <- 2
           stepBounds <- matrix(c(0,Inf),ncol=2,nrow=2*nbStates,byrow=TRUE)
           parNames <- c("mean","sd")
         },
         "weibull"={
           parSize[1] <- 2
           stepBounds <- matrix(c(0,Inf),ncol=2,nrow=2*nbStates,byrow=TRUE)
           parNames <- c("mean","sd")
         },
         "exp"={
           parSize[1] <- 1
           stepBounds <- matrix(c(0,Inf),ncol=2,nrow=nbStates,byrow=TRUE)
           parNames <- c("rate")
         })
  switch(angleDist,
         "vm"={
           if(estAngleMean) {
             parSize[2] <- 2
             angleBounds <- matrix(c(rep(c(-pi,pi),nbStates),rep(c(0,Inf),nbStates)),
                                   ncol=2,byrow=TRUE)
             parNames <- c(parNames,"mean","concentration")
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(c(0,Inf),ncol=2,nrow=nbStates,byrow=TRUE)
             parNames <- c(parNames,"concentration")
           }
         },
         "wrpcauchy"={
           if(estAngleMean) {
             parSize[2] <- 2
             angleBounds <- matrix(c(rep(c(-pi,pi),nbStates),rep(c(0,1),nbStates)),
                                   ncol=2,byrow=TRUE)
             parNames <- c(parNames,"location","concentration")
           }
           else {
             parSize[2] <- 1
             angleBounds <- matrix(c(0,1),ncol=2,nrow=nbStates,byrow=TRUE)
             parNames <- c(parNames,"concentration")
           }
         })

  bounds <- rbind(stepBounds,angleBounds)
  return(list(parSize=parSize,bounds=bounds,parNames=parNames))
}
