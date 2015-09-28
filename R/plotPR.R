
#' Generic plotPR method
#' @param m Fitted model
#' @export
plotPR <- function(m)
  UseMethod("plotPR")

#' Plot pseudo-residuals
#'
#' Plots qq-plots and time series of the pseudo-residuals
#'
#' @method plotPR moveHMM
#'
#' @param m A \code{\link{moveHMM}} object
#'
#' @examples
#' m <- ex$m # moveHMM object (as returned by \code{\link{fitHMM}})
#'
#' plotPR(m)
#'
#' @export

plotPR.moveHMM <- function(m)
{
  cat("Computing pseudo-residuals... ")
  pr <- pseudoRes(m)
  cat("DONE\n")

  par(mfrow=c(2,2))

  # steps qq-plot
  qqStep <- qqnorm(pr$stepRes,plot=FALSE)
  limInf <- min(min(qqStep$x,na.rm=T),min(qqStep$y,na.rm=T))
  limSup <- max(max(qqStep$x,na.rm=T),max(qqStep$y,na.rm=T))
  qqnorm(pr$stepRes,main="Steps pseudo-residuals",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
  abline(0,1,lwd=2)

  # angles qq-plot
  qqAngle <- qqnorm(pr$angleRes,plot=FALSE)
  limInf <- min(min(qqAngle$x,na.rm=T),min(qqAngle$y,na.rm=T))
  limSup <- max(max(qqAngle$x,na.rm=T),max(qqAngle$y,na.rm=T))
  qqnorm(pr$angleRes,main="Angles pseudo-residuals",col="red",xlim=c(limInf,limSup),ylim=c(limInf,limSup))
  abline(0,1,lwd=2)

  # time series
  plot(pr$stepRes,type="l",xlab="Observation index",ylab="Steps pseudo-residuals")
  plot(pr$angleRes,type="l",xlab="Observation index",ylab="Angles pseudo-residuals")

  # back to default
  par(mfrow=c(1,1))
}
