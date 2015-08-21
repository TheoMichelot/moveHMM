
#' Plot moveHMM
#' @method plot moveHMM
#'
#' @param x Object moveHMM
#' @param ask If TRUE, the execution pauses between each plot.
#' @param animals Vector of indices of animals for which information will be plotted.
#' Default : NULL ; all animals are plotted.
#' @param breaks Histogram parameter. See hist documentation.
#' @param hist.ylim Parameter ylim for the step length histograms. See hist documentation.
#' Default : NULL ; the function sets default values.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' m <- ex$mod # moveHMM object, as returned by fitHMM
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)

plot.moveHMM <- function(x,ask=TRUE,animals=NULL,breaks="Sturges",hist.ylim=NULL,...)
{
  m <- x
  nbAnimals <- length(unique(m$data$ID))
  nbStates <- ncol(m$mle$stepPar)

  stepFun <- paste("d",m$stepDist,sep="")
  if(m$angleDist!="none")
    angleFun <- paste("d",m$angleDist,sep="")

  # check arguments
  if(is.null(animals))
    animals <- 1:nbAnimals
  if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0)
    stop("Check animals argument.")

  if(!is.null(hist.ylim) & length(hist.ylim)!=2)
    stop("hist.ylim needs to be a vector of two values (ymin,ymax)")

  # proportion of each state in the states sequence returned by the Viterbi algorithm
  w <- rep(NA,nbStates)
  for(state in 1:nbStates)
    w[state] <- length(which(m$states==state))/length(m$states)

  if(m$conditions$zeroInflation) {
    zeromass <- m$mle$stepPar[nrow(m$mle$stepPar),]
    m$mle$stepPar <- m$mle$stepPar[-nrow(m$mle$stepPar),]
  }

  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=ask)

  if(m$angleDist=="none") { # if step only
    for(zoo in animals) {
      ID <- unique(m$data$ID)[zoo]
      ind <- which(m$data$ID==ID)

      # Histogram of step lengths
      if(is.null(hist.ylim)) { # default
        ymin <- 0
        h <- hist(m$data$step[ind],plot=F)
        ymax <- 1.5*max(h$density)
      }
      else {
        ymin <- hist.ylim[1]
        ymax <- hist.ylim[2]
      }
      hist(m$data$step[ind],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
           col="lightgrey",border="white",breaks=breaks)
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      grid <- seq(0,max(m$data$step[ind],na.rm=T),length=1000)
      for(state in 1:nbStates) {
        # Constitute the lists of state-dependent parameters for the step
        stepArgs <- list(grid)

        for(j in 1:nrow(m$mle$stepPar))
          stepArgs[[j+1]] <- m$mle$stepPar[j,state]

        # conversion between mean/sd and shape/scale if necessary
        if(m$stepDist=="gamma") {
          shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
          scale <- stepArgs[[3]]^2/stepArgs[[2]]
          stepArgs[[2]] <- shape
          stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        # add state-dependent densities to the histogram
        # (weighted by the proportion of each state in the Viterbi states sequence)
        if(m$conditions$zeroInflation)
          lines(grid,(1-zeromass[state])*w[state]*do.call(stepFun,stepArgs),col=state+2,lwd=2)
        else
          lines(grid,w[state]*do.call(stepFun,stepArgs),col=state+2,lwd=2)
      }
    }
  }
  else { # if step + angle
    for(zoo in animals) {
      ID <- unique(m$data$ID)[zoo]
      ind <- which(m$data$ID==ID)
      x <- m$data$x[ind]
      y <- m$data$y[ind]
      states <- m$states[ind]

      # Map of the track, colored by states
      plot(x[1],y[1],xlim=c(min(x,na.rm=T),max(x,na.rm=T)),ylim=c(min(y,na.rm=T),max(y,na.rm=T)),
           pch=18,xlab="x",ylab="y")
      for(i in 2:length(x)) {
        points(x[i],y[i],pch=16,col=states[i-1]+2,cex=0.6)
        segments(x0=x[i-1],y0=y[i-1],x1=x[i],y1=y[i],col=states[i-1]+2,lwd=1.3)
      }
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)

      # Histogram of step lengths
      if(is.null(hist.ylim)) { # default
        ymin <- 0
        h <- hist(m$data$step[ind],plot=F)
        ymax <- 1.5*max(h$density)
      }
      else {
        ymin <- hist.ylim[1]
        ymax <- hist.ylim[2]
      }
      hist(m$data$step[ind],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
           col="lightgrey",border="white",breaks=breaks)
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      grid <- seq(0,max(m$data$step[ind],na.rm=T),length=1000)
      for(state in 1:nbStates) {
        # Constitute the lists of state-dependent parameters for the step
        stepArgs <- list(grid)

        for(j in 1:nrow(m$mle$stepPar))
          stepArgs[[j+1]] <- m$mle$stepPar[j,state]

        # conversion between mean/sd and shape/scale if necessary
        if(m$stepDist=="gamma") {
          shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
          scale <- stepArgs[[3]]^2/stepArgs[[2]]
          stepArgs[[2]] <- shape
          stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }
        # add state-dependent densities to the histogram
        # (weighted by the proportion of each state in the Viterbi states sequence)
        if(m$conditions$zeroInflation)
          lines(grid,(1-zeromass[state])*w[state]*do.call(stepFun,stepArgs),col=state+2,lwd=2)
        else
          lines(grid,w[state]*do.call(stepFun,stepArgs),col=state+2,lwd=2)
      }

      # Histogram of turning angles
      h <- hist(m$data$angle[ind],plot=F) # to determine ylim
      ymax <- 1.5*max(h$density)
      hist(m$data$angle[ind],prob=T,main="",ylim=c(0,ymax),xlab="turning angle (radians)",
           col="lightgrey",border="white",breaks=breaks,xaxt="n")
      axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
           labels = expression(-pi, -pi/2, 0, pi/2, pi))
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      grid <- seq(-pi,pi,length=1000)

      for(state in 1:nbStates) {
        # Constitute the lists of state-dependent parameters for the angle
        angleArgs <- list(grid)

        for(j in 1:nrow(m$mle$anglePar))
          angleArgs[[j+1]] <- m$mle$anglePar[j,state]

        # add state-dependent densities to the histogram
        # (weighted by the proportion of each state in the Viterbi states sequence)
        if(m$conditions$zeroInflation)
          lines(grid,(1-zeromass[state])*w[state]*do.call(angleFun,angleArgs),col=state+2,lwd=2)
        else
          lines(grid,w[state]*do.call(angleFun,angleArgs),col=state+2,lwd=2)
      }
    }
  }

  # plot the transition probabilities as functions of the covariates
  par(mfrow=c(nbStates,nbStates))
  par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right

  covsCol <- which(names(m$data)!="ID" & names(m$data)!="x" & names(m$data)!="y" &
                     names(m$data)!="step" & names(m$data)!="angle")
  allCovs <- m$data[,covsCol]

  if(nrow(m$mle$beta)>1) {
    for(cov in 2:nrow(m$mle$beta)) {
      inf <- min(allCovs[,cov],na.rm=T)
      sup <- max(allCovs[,cov],na.rm=T)

      meanCovs <- colSums(allCovs)/nrow(allCovs)
      desMat <- matrix(rep(meanCovs,100),ncol=length(meanCovs),byrow=TRUE)

      desMat[,cov] <- seq(inf,sup,length=100)

      trMat <- trMatrix_rcpp(nbStates,m$mle$beta,desMat)

      for(i in 1:nbStates)
        for(j in 1:nbStates)
          plot(desMat[,cov],trMat[i,j,],type="l",ylim=c(0,1),xlab=names(allCovs)[cov],
               ylab=paste(i,"->",j))

      mtext("Transition probabilities",side=3,outer=TRUE,padj=2)
    }
  }

  # set the graphical parameters back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
  par(ask=FALSE)
}
