
#' Plot function for a moveData object
#'
#' @param data An object moveData
#' @param compact TRUE for a compact plot (all individuals at once), FALSE otherwise (default -- one
#' individual at a time)

plot.moveData <- function(data,compact=FALSE)
{
  # check arguments
  if(length(data)<1) stop("The data input is empty.")
  if(is.null(data$ID) | is.null(data$x) | is.null(data$y) | is.null(data$step))
    stop("Missing field(s) in data.")

  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=TRUE)

  nbAnimals <- length(unique(data$ID))

  if(is.null(data$angle)) # only step length is provided
  {
    par(mfrow=c(1,2))
    for(zoo in 1:nbAnimals) {
      ID <- unique(data$ID)[zoo]
      step <- data$step[which(data$ID==ID)]
      plot(step,type="l",xlab="t",ylab="step length",
           ylim=c(0,max(step,na.rm=T)))
      hist(step,xlab="step length",main="")
      mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
    }
  }
  else # step length and turning angle are provided
  {
    if(!compact) {
      for(zoo in 1:nbAnimals) {
        ID <- unique(data$ID)[zoo]
        x <- data$x[which(data$ID==ID)]
        y <- data$y[which(data$ID==ID)]
        step <- data$step[which(data$ID==ID)]
        angle <- data$angle[which(data$ID==ID)]

        par(mfrow=c(1,1))
        plot(x,y,type="o",lwd=1.3,xlab="x",ylab="y",pch=20)
        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)

        par(mfrow=c(2,2))
        plot(step,type="l",xlab="t",ylab="step length",
             ylim=c(0,max(step,na.rm=T)))
        plot(angle,type="l",xlab="t",ylab="turning angle",
             ylim=c(-pi,pi))
        abline(h=c(-pi,0,pi),lty=2)

        hist(step,xlab="step length",main="")
        hist(angle,xlab="turning angle",main="")

        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }
    }
    else {
      par(mfrow = c(1,1))
      xmin <- Inf; xmax <- -Inf
      ymin <- Inf; ymax <- -Inf
      stepmax <- -Inf; nbObs <- -Inf
      for(zoo in 1:nbAnimals) {
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

      ID <- unique(data$ID)[1]
      x <- data$x[which(data$ID==ID)]
      y <- data$y[which(data$ID==ID)]
      plot(x,y,type="o",pch=20,lwd=1.3,col=2, cex=0.5,
           xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="x",ylab="y")

      for(zoo in 2:length(data)) {
        ID <- unique(data$ID)[zoo]
        x <- data$x[which(data$ID==ID)]
        y <- data$y[which(data$ID==ID)]
        points(x,y,type="o",pch=20,lwd=1.3,col=zoo+1,cex=0.5)
      }

      par(mfrow=c(2,2))
      for(zoo in 1:nbAnimals) {
        ID <- unique(data$ID)[zoo]
        step <- data$step[which(data$ID==ID)]
        angle <- data$angle[which(data$ID==ID)]

        plot(step,type="l",xlab="t",ylab="step length",
             ylim=c(0,max(step,na.rm=T)))
        plot(angle,type="l",xlab="t",ylab="turning angle",
             ylim=c(-pi,pi))
        abline(h=c(-pi,0,pi),lty=2)

        hist(step,xlab="step length",main="")
        hist(angle,xlab="turning angle",main="")

        mtext(paste("Animal ID :",ID),side=3,outer=TRUE,padj=2)
      }
    }
  }

  par(ask=FALSE)
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # back to default layout
}
