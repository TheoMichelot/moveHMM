
#' Plot moveHMM
#'
#' @param m Object moveHMM
plot.moveHMM <- function(m)
{
  nbAnimals <- length(unique(m$data$ID))

  stepFun <- paste("d",m$stepDist,sep="")
  angleFun <- paste("d",m$angleDist,sep="")

  if(!is.null(m$angleMean))
    m$mle$anglePar <- rbind(m$angleMean,m$mle$anglePar)

  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=TRUE)
  for(zoo in 1:1) {
    ID <- unique(m$data$ID)[zoo]
    ind <- which(m$data$ID==ID)
    x <- m$data$x[ind]
    y <- m$data$y[ind]
    states <- m$states[ind]

    # Map of the track, colored by states
    plot(x[1],y[1],xlim=c(min(x,na.rm=T),max(x,na.rm=T)),ylim=c(min(y,na.rm=T),max(y,na.rm=T)),
         pch=18,xlab="x",ylab="y")
    for(i in 2:length(x)) {
      points(x[i],y[i],pch=16,col=states[i-1]+1,cex=0.6)
      segments(x0=x[i-1],y0=y[i-1],x1=x[i],y1=y[i],col=states[i-1]+1,lwd=1.3)
    }
    mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)

    # Histogram of step lengths
    h <- hist(m$data$step[ind],plot=F) # to choose ylim
    ymax <- 1.5*max(h$density)
    hist(m$data$step[ind],prob=T,main="",ylim=c(0,ymax),xlab="step length")
    mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
    grid <- seq(0,max(m$data$step[ind],na.rm=T),length=1000)
    for(state in 1:nbStates) {
      # Constitute the lists of state-dependent parameters for the step and angle
      stepArgs <- list(grid)

      if(nrow(m$mle$stepPar)==1) stepArgs[[2]] <- m$mle$stepPar[state]
      else {
        for(j in 1:nrow(m$mle$stepPar))
          stepArgs[[j+1]] <- m$mle$stepPar[j,state]
      }
      # conversion between mean/sd and shape/scale if necessary
      if(m$stepDist=="weibull" | m$stepDist=="gamma") {
        shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
        scale <- stepArgs[[3]]^2/stepArgs[[2]]
        stepArgs[[2]] <- shape
        if(m$stepDist=="gamma") stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        else stepArgs[[3]] <- scale # dweibull expects scale
      }
      lines(grid,m$mle$delta[state]*do.call(stepFun,stepArgs),col=state+1,lwd=2)
    }

    # Histogram of turning angles
    if(m$angleDist!="NULL") {
      h <- hist(m$data$angle[ind],plot=F) # to choose ylim
      ymax <- 1.5*max(h$density)
      hist(m$data$angle[ind],prob=T,main="",ylim=c(0,ymax),xlab="turning angle")
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      grid <- seq(-pi,pi,length=1000)

      for(state in 1:nbStates) {
        angleArgs <- list(grid)

        if(nrow(m$mle$anglePar)==1) angleArgs[[2]] <- m$mle$anglePar[state]
        else {
          for(j in 1:nrow(m$mle$anglePar))
            angleArgs[[j+1]] <- m$mle$anglePar[j,state]
        }
        lines(grid,m$mle$delta[state]*do.call(angleFun,angleArgs),col=state+1,lwd=2)
      }
    }
  }
}
