
#' Plot \code{moveData}
#' @method plot moveData
#'
#' @param x An object \code{moveData}
#' @param compact \code{TRUE} for a compact plot (all individuals at once), \code{FALSE} otherwise
#' (default -- one individual at a time).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' data <- ex$data # moveData object, as returned by prepData or simData
#'
#' plot(data,compact=TRUE,breaks=20,ask=FALSE)

plot.moveData <- function(x,compact=FALSE,ask=TRUE,breaks="Sturges",...)
{
  data <- x

  # check arguments
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$ID) | is.null(data$x) | is.null(data$step))
    stop("Missing field(s) in data.")

  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=ask)

  nbAnimals <- length(unique(data$ID))

  if(is.null(data$angle) | length(which(!is.na(data$angle)))==0) # only step length is provided
  {
    par(mfrow=c(1,2))
    for(zoo in 1:nbAnimals) {
      ID <- unique(data$ID)[zoo]
      step <- data$step[which(data$ID==ID)]
      # step length time series
      plot(step,type="l",xlab="t",ylab="step length",
           ylim=c(0,max(step,na.rm=T)))
      # step length histogram
      hist(step,xlab="step length",main="",col="lightblue",border="white",breaks=breaks)
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
    }
  }
  else # step length and turning angle are provided
  {
    if(!compact) { # tracks are plotted on a separate map for each animal
      for(zoo in 1:nbAnimals) {
        ID <- unique(data$ID)[zoo]
        x <- data$x[which(data$ID==ID)]
        y <- data$y[which(data$ID==ID)]
        step <- data$step[which(data$ID==ID)]
        angle <- data$angle[which(data$ID==ID)]

        par(mfrow=c(1,1))
        # map of the animal's track
        plot(x,y,type="o",lwd=1.3,xlab="x",ylab="y",pch=20)
        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)

        # step and angle time series
        par(mfrow=c(2,2))
        plot(step,type="l",xlab="t",ylab="step length",
             ylim=c(0,max(step,na.rm=T)))
        plot(angle,type="l",xlab="t",ylab="turning angle (radians)",
             ylim=c(-pi,pi),yaxt="n")
        axis(2, at = c(-pi, -pi/2, 0, pi/2, pi),
             labels = expression(-pi, -pi/2, 0, pi/2, pi))
        abline(h=c(-pi,0,pi),lty=2)

        # step and angle histograms
        hist(step,xlab="step length",main="", col="lightgrey",border="white",breaks=breaks)
        hist(angle,xlab="turning angle (radians)",main="", col="lightgrey",border="white",
             breaks=breaks,xaxt="n")
        axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
             labels = expression(-pi, -pi/2, 0, pi/2, pi))

        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }
    }
    else { # tracks are plotted on a single map for all animals
      par(mfrow = c(1,1))
      xmin <- Inf; xmax <- -Inf
      ymin <- Inf; ymax <- -Inf
      stepmax <- -Inf; nbObs <- -Inf
      for(zoo in 1:nbAnimals) {
        # compute the x and y bounds to draw the map
        ID <- unique(data$ID)[zoo]
        x <- data$x[which(data$ID==ID)]
        y <- data$y[which(data$ID==ID)]
        step <- data$step[which(data$ID==ID)]

        if(min(x,na.rm=T)<xmin) xmin <- min(x,na.rm=T) # na.rm=T to ignore the NAs
        if(min(y,na.rm=T)<ymin) ymin <- min(y,na.rm=T)
        if(max(x,na.rm=T)>xmax) xmax <- max(x,na.rm=T)
        if(max(y,na.rm=T)>ymax) ymax <- max(y,na.rm=T)
        if(max(step,na.rm=T)>stepmax) stepmax <- max(step,na.rm=T)
        if(length(x)>nbObs) nbObs <- length(x)
      }

      if(nbAnimals>6)
        colors <- rainbow(nbAnimals) # to make sure that all colors are distinct
      else
        colors <- c(2,3,4,5,6,7)

      ID <- unique(data$ID)[1]
      x <- data$x[which(data$ID==ID)]
      y <- data$y[which(data$ID==ID)]
      # plot the first animal's track
      plot(x,y,type="o",pch=20,lwd=1.3,col=colors[1], cex=0.5,
           xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="x",ylab="y")

      # add each other animal's track to the map
      for(zoo in 2:nbAnimals) {
        ID <- unique(data$ID)[zoo]
        x <- data$x[which(data$ID==ID)]
        y <- data$y[which(data$ID==ID)]
        points(x,y,type="o",pch=20,lwd=1.3,col=colors[zoo],cex=0.5)
      }

      par(mfrow=c(2,2))
      for(zoo in 1:nbAnimals) {
        ID <- unique(data$ID)[zoo]
        step <- data$step[which(data$ID==ID)]
        angle <- data$angle[which(data$ID==ID)]

        # step and angle time series
        plot(step,type="l",xlab="t",ylab="step length",
             ylim=c(0,max(step,na.rm=T)))
        plot(angle,type="l",xlab="t",ylab="turning angle (radians)",
             ylim=c(-pi,pi),yaxt="n")
        axis(2, at = c(-pi, -pi/2, 0, pi/2, pi),
             labels = expression(-pi, -pi/2, 0, pi/2, pi))
        abline(h=c(-pi,0,pi),lty=2)

        # step and angle histograms
        hist(step,xlab="step length",main="", col="lightgrey",border="white",breaks=breaks)
        hist(angle,xlab="turning angle (radians)",main="", col="lightgrey",border="white",breaks=breaks,
             xaxt="n")
        axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
             labels = expression(-pi, -pi/2, 0, pi/2, pi))

        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }
    }
  }

  # set graphical parameters back to default
  par(ask=FALSE)
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2))
}
