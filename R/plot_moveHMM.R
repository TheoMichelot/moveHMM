
#' Plot \code{moveHMM}
#' @method plot moveHMM
#'
#' @param x Object \code{moveHMM}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default : \code{NULL} ; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim Parameter \code{ylim} for the step length histograms. See \code{hist} documentation.
#' Default : \code{NULL} ; the function sets default values.
#' @param compactHist If \code{TRUE}, the function only plots histograms of all observations for steps
#' and angles, with the fitted densities (no map, and no individual-specific plot). Default : \code{FALSE}.
#' @param sepStates If TRUE, the data is split by states in the histograms. Default : \code{FALSE}.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' m <- ex$m # moveHMM object, as returned by fitHMM
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)

plot.moveHMM <- function(x,animals=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,compactHist=FALSE,
                         sepStates=FALSE,...)
{
  m <- x
  nbAnimals <- length(unique(m$data$ID))
  nbStates <- ncol(m$mle$stepPar)

  stepFun <- paste("d",m$stepDist,sep="")
  if(m$angleDist!="none")
    angleFun <- paste("d",m$angleDist,sep="")

  # define animals to be plotted
  if(is.null(animals)) # all animals are plotted
    animalsInd <- 1:nbAnimals
  else {
    if(is.character(animals)) { # animals' IDs provided
      animalsInd <- NULL
      for(zoo in 1:length(animals)) {
        if(length(which(unique(m$data$ID)==animals[zoo]))==0) # ID not found
          stop("Check animals argument.")

        animalsInd <- c(animalsInd,which(unique(m$data$ID)==animals[zoo]))
      }
    }

    if(is.numeric(animals)) { # animals' indices provided
      if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
        stop("Check animals argument.")

      animalsInd <- animals
    }
  }

  # if compactHist, only one histogram is plotted
  if(compactHist)
    animalsInd <- 1

  # check arguments
  if(!is.null(hist.ylim) & length(hist.ylim)!=2)
    stop("hist.ylim needs to be a vector of two values (ymin,ymax)")

  # compute most probable states sequence using Viterbi
  if(nbStates>1) {
    cat("Decoding states sequence... ")
    vitStates <- viterbi(m)
    cat("DONE\n")
  } else
    vitStates <- rep(0,nrow(m$data))

  if(sepStates | nbStates==1)
    w <- rep(1,nbStates)
  else {
    # proportion of each state in the states sequence returned by the Viterbi algorithm
    w <- rep(NA,nbStates)
    for(state in 1:nbStates)
      w[state] <- length(which(vitStates==state))/length(vitStates)
  }

  if(m$conditions$zeroInflation) {
    zeromass <- m$mle$stepPar[nrow(m$mle$stepPar),]
    m$mle$stepPar <- m$mle$stepPar[-nrow(m$mle$stepPar),]
  }

  # text for legends
  legText <- NULL
  for(i in 1:nbStates)
    legText <- c(legText,paste("State",i))

  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=ask)

  if(m$angleDist=="none") { # if step only
    for(zoo in animalsInd) {
      if(compactHist) {
        ID <- "all animals"
        ind <- 1:nrow(m$data)
      } else {
        ID <- unique(m$data$ID)[zoo]
        ind <- which(m$data$ID==ID)
      }

      # Histogram of step lengths
      if(!sepStates) {
        if(is.null(hist.ylim)) { # default
          ymin <- 0
          h <- hist(m$data$step[ind],plot=F,breaks=breaks)
          ymax <- 1.3*max(h$density)
        } else {
          ymin <- hist.ylim[1]
          ymax <- hist.ylim[2]
        }

        hist(m$data$step[ind],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
             col="lightgrey",border="white",breaks=breaks)

        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }

      grid <- seq(0,max(m$data$step[ind],na.rm=T),length=1000)
      for(state in 1:nbStates) {
        if(sepStates) {
          stateInd <- ind[which(vitStates[ind]==state)]

          if(is.null(hist.ylim)) { # default
            ymin <- 0
            h <- hist(m$data$step[stateInd],plot=F,breaks=breaks)
            ymax <- 1.3*max(h$density)
          }
          else {
            ymin <- hist.ylim[1]
            ymax <- hist.ylim[2]
          }

          hist(m$data$step[stateInd],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
               col="lightgrey",border="white",breaks=breaks)

          mtext(paste("Animal ID :",ID," -- State",state),side=3,outer=TRUE,padj=2)
        }

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
          lines(grid,(1-zeromass[state])*w[state]*do.call(stepFun,stepArgs),col=state+1,lwd=2)
        else
          lines(grid,w[state]*do.call(stepFun,stepArgs),col=state+1,lwd=2)
      }

      # add a legend
      if(!sepStates)
        legend("top",legText,lwd=rep(2,nbStates),col=c(2:(nbStates+1)),bty="n")
    }
  }
  else { # if step + angle
    for(zoo in animalsInd) {
      if(compactHist) {
        ID <- "all animals"
        ind <- 1:nrow(m$data)
      } else {
        ID <- unique(m$data$ID)[zoo]
        ind <- which(m$data$ID==ID)
        x <- m$data$x[ind]
        y <- m$data$y[ind]
        states <- vitStates[ind]

        # Map of the track, colored by states
        plot(x[1],y[1],xlim=c(min(x,na.rm=T),max(x,na.rm=T)),ylim=c(min(y,na.rm=T),max(y,na.rm=T)),
             pch=18,xlab="x",ylab="y")
        for(i in 2:length(x)) {
          points(x[i],y[i],pch=16,col=states[i-1]+1,cex=0.6)
          segments(x0=x[i-1],y0=y[i-1],x1=x[i],y1=y[i],col=states[i-1]+1,lwd=1.3)
        }
        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }

      # Histogram of step lengths
      if(!sepStates) {
        if(is.null(hist.ylim)) { # default
          ymin <- 0
          h <- hist(m$data$step[ind],plot=F,breaks=breaks)
          ymax <- 1.3*max(h$density)
        }
        else {
          ymin <- hist.ylim[1]
          ymax <- hist.ylim[2]
        }

        hist(m$data$step[ind],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
             col="lightgrey",border="white",breaks=breaks)
        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }

      grid <- seq(0,max(m$data$step[ind],na.rm=T),length=1000)
      for(state in 1:nbStates) {
        if(sepStates) {
          stateInd <- ind[which(vitStates[ind]==state)]

          if(is.null(hist.ylim)) { # default
            ymin <- 0
            h <- hist(m$data$step[stateInd],plot=F,breaks=breaks)
            ymax <- 1.3*max(h$density)
          }
          else {
            ymin <- hist.ylim[1]
            ymax <- hist.ylim[2]
          }

          hist(m$data$step[stateInd],prob=T,main="",ylim=c(ymin,ymax),xlab="step length",
               col="lightgrey",border="white",breaks=breaks)

          mtext(paste("Animal ID :",ID," -- State",state),side=3,outer=TRUE,padj=2)
        }

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
          lines(grid,(1-zeromass[state])*w[state]*do.call(stepFun,stepArgs),col=state+1,lwd=2)
        else
          lines(grid,w[state]*do.call(stepFun,stepArgs),col=state+1,lwd=2)
      }

      # add a legend
      if(!sepStates)
        legend("top",legText,lwd=rep(2,nbStates),col=c(2:(nbStates+1)),bty="n")

      # Histogram of turning angles
      if(!sepStates) {
        h <- hist(m$data$angle[ind],plot=F,breaks=breaks) # to determine ylim
        ymax <- 1.3*max(h$density)

        hist(m$data$angle[ind],prob=T,main="",ylim=c(0,ymax),xlab="turning angle (radians)",
             col="lightgrey",border="white",breaks=seq(-pi,pi,length=length(h$breaks)),xaxt="n")
        axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
             labels = expression(-pi, -pi/2, 0, pi/2, pi))
        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }

      grid <- seq(-pi,pi,length=1000)

      for(state in 1:nbStates) {
        if(sepStates) {
          stateInd <- ind[which(vitStates[ind]==state)]

          h <- hist(m$data$angle[stateInd],plot=F,breaks=breaks) # to determine ylim
          ymax <- 1.3*max(h$density)

          hist(m$data$angle[stateInd],prob=T,main="",ylim=c(0,ymax),xlab="turning angle (radians)",
               col="lightgrey",border="white",breaks=seq(-pi,pi,length=length(h$breaks)),xaxt="n")
          axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
               labels = expression(-pi, -pi/2, 0, pi/2, pi))
          mtext(paste("Animal ID :",ID," -- State",state),side=3,outer=TRUE,padj=2)
        }

        # Constitute the lists of state-dependent parameters for the angle
        angleArgs <- list(grid)

        for(j in 1:nrow(m$mle$anglePar))
          angleArgs[[j+1]] <- m$mle$anglePar[j,state]

        # add state-dependent densities to the histogram
        # (weighted by the proportion of each state in the Viterbi states sequence)
        if(m$conditions$zeroInflation)
          lines(grid,(1-zeromass[state])*w[state]*do.call(angleFun,angleArgs),col=state+1,lwd=2)
        else
          lines(grid,w[state]*do.call(angleFun,angleArgs),col=state+1,lwd=2)
      }

      # add a legend
      if(!sepStates)
        legend("top",legText,lwd=rep(2,nbStates),col=c(2:(nbStates+1)),bty="n")
    }
  }

  # plot the transition probabilities as functions of the covariates
  if(!compactHist & nbStates>1) {
    par(mfrow=c(nbStates,nbStates))
    par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right

    rawCovs <- m$rawCovs
    gridLength <- 100

    if(nrow(m$mle$beta)>1) {
      for(cov in 1:ncol(m$rawCovs)) {
        inf <- min(rawCovs[,cov],na.rm=T)
        sup <- max(rawCovs[,cov],na.rm=T)

        # mean values of each covariate
        meanCovs <- colSums(rawCovs)/nrow(rawCovs)

        # set all covariates to their mean, except for "cov"
        # (which takes a grid of values from inf to sup)
        tempCovs <- data.frame(rep(meanCovs[1],gridLength))
        if(length(meanCovs)>1)
           for(i in 2:length(meanCovs))
             tempCovs <- cbind(tempCovs,rep(meanCovs[i],gridLength))

        tempCovs[,cov] <- seq(inf,sup,length=gridLength)
        colnames(tempCovs) <- colnames(rawCovs)

        desMat <- model.matrix(m$conditions$formula,data=tempCovs)

        # check that the current covariate (cov) is included in the model
        used <- FALSE
        for(i in 2:ncol(desMat)) {
          c <- desMat[,i]
          if(length(which(c!=mean(c)))>0)
            used <- TRUE
        }

        if(used) {
          trMat <- trMatrix_rcpp(nbStates,m$mle$beta,desMat)

          for(i in 1:nbStates)
            for(j in 1:nbStates)
              plot(tempCovs[,cov],trMat[i,j,],type="l",ylim=c(0,1),xlab=names(rawCovs)[cov],
                   ylab=paste(i,"->",j))

          mtext("Transition probabilities",side=3,outer=TRUE,padj=2)
        }
      }
    }
  }

  # set the graphical parameters back to default
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # bottom, left, top, right
  par(ask=FALSE)
}
