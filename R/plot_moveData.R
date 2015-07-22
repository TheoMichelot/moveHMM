
#' Plot function for a moveData object
#'
#' @param data An object moveData
#' @param compact TRUE for a compact plot (all individuals at once), FALSE otherwise (default -- one
#' individual at a time)

plot.moveData <- function(data,compact=FALSE)
{
  par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
  par(ask=T)
  if(!compact) {
    for(i in 1:length(data)) {
      par(mfrow=c(1,1))
      plot(data[[i]]$x,data[[i]]$y,type="o",lwd=1.3,xlab="x",ylab="y",pch=20)
      mtext(paste("Animal ID :",data[[i]]$ID),side=3,outer=TRUE,padj=2)

      par(mfrow=c(2,2))
      plot(data[[i]]$step,type="l",xlab="t",ylab="step length",
           ylim=c(0,max(data[[i]]$step,na.rm=T)))
      plot(data[[i]]$angle,type="l",xlab="t",ylab="turning angle",
           ylim=c(-pi,pi))
      abline(h=c(-pi,0,pi),lty=2)

      hist(data[[i]]$step,xlab="step length",main="")
      hist(data[[i]]$angle,xlab="turning angle",main="")

      mtext(paste("Animal ID :",data[[i]]$ID),side=3,outer=TRUE,padj=2)
    }
  }
  else {
    par(mfrow = c(1,1))
    xmin <- Inf; xmax <- -Inf
    ymin <- Inf; ymax <- -Inf
    stepmax <- -Inf; nbObs <- -Inf
    for(i in 1:length(data)) {
      if(min(data[[i]]$x,na.rm=T)<xmin) xmin <- min(data[[i]]$x,na.rm=T) # na.rm=T to ignore the NAs
      if(min(data[[i]]$y,na.rm=T)<ymin) ymin <- min(data[[i]]$y,na.rm=T)
      if(max(data[[i]]$x,na.rm=T)>xmax) xmax <- max(data[[i]]$x,na.rm=T)
      if(max(data[[i]]$y,na.rm=T)>ymax) ymax <- max(data[[i]]$y,na.rm=T)
      if(max(data[[i]]$step,na.rm=T)>stepmax) stepmax <- max(data[[i]]$step,na.rm=T)
      if(length(data[[i]]$x)>nbObs) nbObs <- length(data[[i]]$x)
    }
    plot(data[[1]]$x,data[[1]]$y,type="o",pch=20,lwd=1.3,col=1, cex=0.5,
         xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="x",ylab="y")
    for(i in 2:length(data)) {
      points(data[[i]]$x,data[[i]]$y,type="o",pch=20,lwd=1.3,col=i,cex=0.5)
    }

    par(mfrow=c(2,2))
    for(i in 1:length(data)) {
      plot(data[[i]]$step,type="l",xlab="t",ylab="step length",
           ylim=c(0,max(data[[i]]$step,na.rm=T)))
      plot(data[[i]]$angle,type="l",xlab="t",ylab="turning angle",
           ylim=c(-pi,pi))
      abline(h=c(-pi,0,pi),lty=2)

      hist(data[[i]]$step,xlab="step length",main="")
      hist(data[[i]]$angle,xlab="turning angle",main="")

      mtext(paste("Animal ID :",data[[i]]$ID),side=3,outer=TRUE,padj=2)
    }

  }
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)) # back to default layout
}
