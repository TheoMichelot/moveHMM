
#' Plot \code{moveData}
#' @method plot moveData
#'
#' @param x An object \code{moveData}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL} ; all animals are plotted.
#' @param compact \code{TRUE} for a compact plot (all individuals at once), \code{FALSE} otherwise
#' (default -- one individual at a time).
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # data is a moveData object (as returned by prepData), automatically loaded with the package
#' data <- example$data
#'
#' plot(data,compact=TRUE,breaks=20,ask=FALSE)
#'
#' @export
#' @importFrom graphics abline axis hist mtext par plot points
#' @importFrom grDevices hcl

plot.moveData <- function(x,animals=NULL,compact=FALSE,ask=TRUE,breaks="Sturges",...)
{
    data <- x
    # check arguments
    if(length(data)<1) stop("The data input is empty.")
    if(is.null(data$ID) | is.null(data$x) | is.null(data$step))
        stop("Missing field(s) in data.")

    nbAnimals <- length(unique(data$ID))

    if(all(data$y==0) & compact) {
        warning("One-dimensional data cannot plotted with the option 'compact'.")
        compact <- FALSE
    }

    ##################################
    ## Define animals to be plotted ##
    ##################################
    if(is.character(animals)) {
        if(any(!animals %in% unique(data$ID))) {
            stop("Check 'animals' argument, ID not found")
        }
        animalsInd <- which(unique(data$ID) %in% animals)
    } else if(is.numeric(animals)) {
        if(min(animals) < 1 | max(animals) > nbAnimals) {
            stop("Check 'animals' argument, index out of bounds")
        }
        animalsInd <- animals
    } else {
        animalsInd <- 1:nbAnimals
    }

    # graphical parameters
    par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
    par(ask=ask)

    if(is.null(data$angle) | all(data$y==0)) {
        # only step length is provided
        for(zoo in animalsInd) {
            ################
            ## Plot track ##
            ################*
            ID <- unique(data$ID)[zoo]
            x <- data$x[which(data$ID==ID)]
            plot(x, type="o",lwd=1.3, xlab="time", ylab="x", pch=20)

            ##########################################
            ## Plot steps time series and histogram ##
            ##########################################
            par(mfrow=c(1,2))
            step <- data$step[which(data$ID==ID)]
            # step length time series
            plot(step,type="l",xlab="t",ylab="step length",
                 ylim=c(0,max(step,na.rm=T)))
            # step length histogram
            hist(step,xlab="step length",main="",col="grey",border="white",breaks=breaks)
            mtext(paste("Animal ID:",ID),side=3,outer=TRUE,padj=2)
        }
    } else {
        # step length and turning angle are provided
        if(compact) {
            ################################
            ## Map of all animals' tracks ##
            ################################
            par(mfrow = c(1,1))
            colors <- getPalette(nbStates = length(animalsInd))

            # determine bounds
            ind <- which(data$ID %in% unique(data$ID)[animalsInd])
            xlim <- range(data$x[ind], na.rm = TRUE)
            ylim <- range(data$y[ind], na.rm = TRUE)

            # plot tracks
            plot(NA, xlim = xlim, ylim = ylim,
                 xlab = "x", ylab = "y", asp = 1)
            for(zoo in animalsInd) {
                ID <- unique(data$ID)[zoo]
                x <- data$x[which(data$ID == ID)]
                y <- data$y[which(data$ID == ID)]
                points(x, y, type = "o", pch = 20, lwd = 1.3,
                       col = colors[zoo], cex = 0.5)
            }
        }

        for(zoo in animalsInd) {
            ID <- unique(data$ID)[zoo]
            step <- data$step[which(data$ID==ID)]
            angle <- data$angle[which(data$ID==ID)]

            if(!compact) {
                # map of the animal's track
                par(mfrow = c(1, 1))
                x <- data$x[which(data$ID==ID)]
                y <- data$y[which(data$ID==ID)]
                plot(x, y, type = "o", lwd = 1.3, xlab = "x", ylab = "y",
                     pch = 20, cex = 0.5, asp = 1)
                mtext(paste("Animal ID:", ID), side = 3, outer = TRUE, padj = 2)
            }

            ##################################
            ## Steps and angles time series ##
            ##################################
            par(mfrow = c(2, 2))
            plot(step, type = "h", xlab = "time", ylab = "step length",
                 ylim=c(0, max(step, na.rm = T)))
            plot(angle, type = "h", xlab = "time", ylab = "turning angle (radians)",
                 ylim = c(-pi, pi), yaxt = "n")
            axis(2, at = c(-pi, -pi/2, 0, pi/2, pi),
                 labels = expression(-pi, -pi/2, 0, pi/2, pi))
            abline(h = c(-pi, 0, pi), lty = 2)

            #################################
            ## Steps and angles histograms ##
            #################################
            hist(step, xlab = "step length", main = "", breaks = breaks,
                 col = "grey", border = "white")

            h <- hist(angle, breaks = breaks, plot = FALSE) # to define the breaks
            hist(angle, xlab = "turning angle (radians)", main = "",
                 col = "grey", border = "white", xaxt="n",
                 breaks = seq(-pi, pi, length = length(h$breaks)))
            axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
                 labels = expression(-pi, -pi/2, 0, pi/2, pi))

            mtext(paste("Animal ID:", ID), side = 3, outer = TRUE, padj = 2)
        }
    }

    # set graphical parameters back to default
    par(ask = FALSE)
    par(mfrow = c(1, 1))
    par(mar = c(5, 4, 4, 2) + 0.1)
}
