
#' Plot \code{moveHMM}
#'
#' Plot the fitted step and angle densities over histograms of the data, transition probabilities
#' as functions of the covariates, and maps of the animals' tracks colored by the decoded states.
#'
#' @method plot moveHMM
#'
#' @param x Object \code{moveHMM}
#' @param animals Vector of indices or IDs of animals for which information will be plotted.
#' Default: \code{NULL}; all animals are plotted.
#' @param ask If \code{TRUE}, the execution pauses between each plot.
#' @param breaks Histogram parameter. See \code{hist} documentation.
#' @param hist.ylim Parameter \code{ylim} for the step length histograms.
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values.
#' @param sepAnimals If \code{TRUE}, the data is split by individuals in the histograms.
#' Default: \code{FALSE}.
#' @param sepStates If \code{TRUE}, the data is split by states in the histograms.
#' Default: \code{FALSE}.
#' @param col Vector or colors for the states (one color per state).
#' @param cumul If \code{TRUE}, the sum of weighted densities is plotted (default).
#' @param plotTracks If \code{TRUE}, the Viterbi-decoded tracks are plotted (default).
#' @param plotCI If \code{TRUE}, confidence intervals are plotted on the transition
#' probabilities (default: FALSE).
#' @param alpha Significance level of the confidence intervals if plotCI=TRUE.
#' Default: 0.95 (i.e. 95\% CIs).
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @details The state-dependent densities are weighted by the frequency of each state in the most
#' probable state sequence (decoded with the function \code{\link{viterbi}}). For example, if the
#' most probable state sequence indicates that one third of observations correspond to the first
#' state, and two thirds to the second state, the plots of the densities in the first state are
#' weighted by a factor 1/3, and in the second state by a factor 2/3.
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plot(m,ask=TRUE,animals=1,breaks=20)
#'
#'
#' @export
#' @importFrom graphics legend lines segments arrows
#' @importFrom grDevices gray
#' @importFrom stats plogis qlogis
#' @importFrom numDeriv grad

plot.moveHMM <- function(x,animals=NULL,ask=TRUE,breaks="Sturges",hist.ylim=NULL,sepAnimals=FALSE,
                         sepStates=FALSE,col=NULL,cumul=TRUE,plotTracks=TRUE,plotCI=FALSE,alpha=0.95,...)
{
    m <- x # the name "x" is for compatibility with the generic method
    nbAnimals <- length(unique(m$data$ID))
    nbStates <- ncol(m$mle$stepPar)

    stepFun <- paste("d",m$conditions$stepDist,sep="")
    if(m$conditions$angleDist!="none")
        angleFun <- paste("d",m$conditions$angleDist,sep="")

    if(!is.null(hist.ylim) & length(hist.ylim)!=2)
        stop("hist.ylim needs to be a vector of two values (ymin,ymax)")

    # prepare colors for the states (used in the maps and for the densities)
    if(!is.null(col) & length(col)!=nbStates) {
        warning("Length of 'col' should be equal to number of states - argument ignored")
        col <- 2:(nbStates+1)
    }
    if(is.null(col) & nbStates<8) {
        # color-blind friendly palette
        pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        col <- pal[1:nbStates]
    }
    if(is.null(col) & nbStates>=8) {
        # to make sure that all colours are distinct (emulate ggplot default palette)
        hues <- seq(15, 375, length = nbStates + 1)
        col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
    }

    if(sepStates | nbStates<2)
        cumul <- FALSE

    ######################
    ## Prepare the data ##
    ######################
    # determine indices of animals to be plotted
    if(is.null(animals)) # all animals are plotted
        animalsInd <- 1:nbAnimals
    else {
        if(is.character(animals)) { # animals' IDs provided
            animalsInd <- NULL
            for(zoo in 1:length(animals)) {
                if(length(which(unique(m$data$ID)==animals[zoo]))==0) # ID not found
                    stop("Check 'animals' argument, ID not found")

                animalsInd <- c(animalsInd,which(unique(m$data$ID)==animals[zoo]))
            }
        }

        if(is.numeric(animals)) { # animals' indices provided
            if(length(which(animals<1))>0 | length(which(animals>nbAnimals))>0) # index out of bounds
                stop("Check 'animals' argument, index out of bounds")

            animalsInd <- animals
        }
    }

    nbAnimals <- length(animalsInd)
    ID <- unique(m$data$ID)[animalsInd]

    # split data by animals if necessary
    if(sepAnimals) {
        stepData <- list()
        angleData <- list()
        for(zoo in 1:nbAnimals) {
            ind <- which(m$data$ID==ID[zoo])
            stepData[[zoo]] <- m$data$step[ind]
            angleData[[zoo]] <- m$data$angle[ind]
        }
    } else {
        ind <- which(m$data$ID %in% ID)
        stepData <- m$data$step[ind]
        angleData <- m$data$angle[ind]
    }

    if(m$conditions$angleDist=="none")
        angleData <- NULL

    x <- list()
    y <- list()
    for(zoo in 1:nbAnimals) {
        ind <- which(m$data$ID==ID[zoo])
        x[[zoo]] <- m$data$x[ind]
        y[[zoo]] <- m$data$y[ind]
    }

    ##################################
    ## States decoding with Viterbi ##
    ##################################
    if(nbStates>1) {
        cat("Decoding states sequence... ")
        states <- viterbi(m)
        cat("DONE\n")
    } else
        states <- rep(1,nrow(m$data))

    ########################################
    ## Plot state-dependent distributions ##
    ########################################
    par(mar=c(5, 4, 4, 2) - c(0, 0, 2, 1)) # bottom, left, top, right
    par(ask = ask)

    distData <- getPlotData(mod, type = "dist")

    # setup line options
    legText <- c(paste("state", 1:nbStates), "total")
    lty <- c(rep(1, nbStates), 2)
    lwd <- c(rep(1, nbStates), 2)
    col <- c(pal[1:nbStates], "black")

    # define ymax for step histogram
    h <- hist(m$data$step, plot = FALSE, breaks = breaks)
    ymax <- 1.3 * max(h$density)
    maxdens <- max(distData$step$total)
    if(maxdens > ymax & maxdens < 1.5 * ymax) {
        ymax <- maxdens
    }

    # step length histogram
    hist(m$data$step, ylim = c(0, ymax), prob = TRUE, main = "",
         xlab = "step length", col = "lightgrey", border = "white",
         breaks = breaks)
    for(i in 1:(nbStates + 1)) {
        lines(distData$step$step, distData$step[,i+1], col = col[i],
              lty = lty[i], lwd = lwd[i])
    }
    legend("top", legText, lwd = lwd, col = col, lty = lty, bty = "n")

    # define ymax and breaks for angle histogram
    h1 <- hist(m$data$angle, plot = FALSE, breaks = breaks)
    breaks <- seq(-pi, pi, length = length(h1$breaks))
    h2 <- hist(m$data$angle, plot = FALSE, breaks = breaks)
    ymax <- 1.3 * max(h2$density)

    # turning angle histogram
    hist(m$data$angle, ylim = c(0, ymax), prob = TRUE, main = "",
         xlab = "turning angle", col = "lightgrey", border = "white",
         breaks = breaks, xaxt = "n")
    axis(1, at = c(-pi, -pi/2, 0, pi/2, pi),
         labels = expression(-pi, -pi/2, 0, pi/2, pi))
    for(i in 1:(nbStates + 1)) {
        lines(distData$angle$angle, distData$angle[,i+1], col = col[i],
              lty = lty[i], lwd = lwd[i])
    }
    legend("top", legText, lwd = lwd, col = col, lty = lty, bty = "n")

    ##################################################
    ## Plot the t.p. as functions of the covariates ##
    ##################################################
    if(nbStates>1) {
        beta <- m$mle$beta

        if(nrow(beta)>1) {
            par(mfrow = c(nbStates, nbStates))
            par(mar = c(5, 4, 4, 2) - c(0, 0, 1.5, 1)) # bottom, left, top, right

            trProbs <- getPlotData(m, type = "tpm", format = "wide")

            # loop over covariates
            par(mfrow=c(nbStates,nbStates))
            par(mar=c(5,4,4,2)-c(0,0,1.5,1)) # bottom, left, top, right
            for(cov in 1:ncol(m$rawCovs)) {
                trProbsCov <- trProbs[[cov]]
                # loop over entries of the transition probability matrix
                for(i in 1:nbStates) {
                    for(j in 1:nbStates) {
                        trName <- paste0("S", i, "toS", j)
                        plot(trProbsCov[,cov], trProbsCov[,trName], type = "l",
                             ylim = c(0, 1), xlab = names(trProbs)[cov],
                             ylab = paste(i, "->", j))

                        # derive confidence intervals using the delta method
                        if(plotCI) {
                            options(warn = -1) # to muffle "zero-length arrow..." warning
                            # plot the confidence intervals
                            arrows(trProbsCov[,cov], trProbsCov[,paste0(trName, ".lci")],
                                   trProbsCov[,cov], trProbsCov[,paste0(trName, ".uci")],
                                   length = 0.025, angle = 90, code = 3,
                                   col = gray(0.5), lwd = 0.7)
                            options(warn = 1)
                        }
                    }
                }

                mtext("Transition probabilities", side = 3, outer = TRUE, padj = 2)
            }
        }
    }

    #################################
    ## Plot maps colored by states ##
    #################################
    if(nbStates>1 & plotTracks) { # no need to plot the map if only one state
        par(mfrow=c(1,1))
        par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right

        for(zoo in 1:nbAnimals) {
            # states for animal 'zoo'
            subStates <- states[which(m$data$ID==ID[zoo])]

            if(!all(y[[zoo]]==0)) { # if 2D data

                plot(x[[zoo]],y[[zoo]],pch=16,col=col[subStates],cex=0.5,asp=1,
                     xlab="x",ylab="y")
                segments(x0=x[[zoo]][-length(x[[zoo]])], y0=y[[zoo]][-length(y[[zoo]])],
                         x1=x[[zoo]][-1], y1=y[[zoo]][-1],
                         col=col[subStates[-length(subStates)]],lwd=1.3)

            } else { # if 1D data

                ymin <- min(x[[zoo]],na.rm=T)
                ymax <- max(x[[zoo]],na.rm=T)

                # first point
                plot(x[[zoo]][1], lwd=1.3, xlim=c(1,length(x[[zoo]])), ylim=c(ymin,ymax),
                     xlab="time", ylab="x", pch=18, col=col[subStates[1]])

                # trajectory
                for(i in 2:length(x[[zoo]])) {
                    points(i,x[[zoo]][i],pch=16,col=col[subStates[i-1]],cex=0.5)
                    segments(x0=i-1,y0=x[[zoo]][i-1],x1=i,y1=x[[zoo]][i],
                             col=col[subStates[i-1]],lwd=1.3)
                }
            }

            mtext(paste("Animal ID:",ID[zoo]),side=3,outer=TRUE,padj=2)
        }
    }

    # set the graphical parameters back to default
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1) # bottom, left, top, right
    par(ask=FALSE)
}
