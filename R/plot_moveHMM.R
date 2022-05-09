
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
#' See \code{hist} documentation. Default: \code{NULL} ; the function sets default values.
#' @param col Vector or colors for the states (one color per state).
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
#' @export
#' @importFrom graphics legend lines segments arrows
#' @importFrom grDevices gray
#' @importFrom stats plogis qlogis
#' @importFrom numDeriv grad

plot.moveHMM <- function(x, animals = NULL, ask = TRUE, breaks = "Sturges", col = NULL,
                         plotTracks = TRUE, plotCI = FALSE, alpha = 0.95, ...) {
    m <- x # the name "x" is for compatibility with the generic method
    nbStates <- ncol(m$mle$stepPar)

    # prepare colours for the states (used in the maps and for the densities)
    if(is.null(col) | (!is.null(col) & length(col) != nbStates)) {
        col <- getPalette(nbStates = nbStates)
    }

    #################################
    ## State decoding with Viterbi ##
    #################################
    if(nbStates > 1) {
        cat("Decoding states sequence... ")
        states <- viterbi(m)
        cat("DONE\n")
    } else {
        states <- rep(1,nrow(m$data))
    }

    ########################################
    ## Plot state-dependent distributions ##
    ########################################
    par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1)) # bottom, left, top, right
    par(ask = ask)

    distData <- getPlotData(m = m, type = "dist")

    # setup line options
    legText <- c(paste("state", 1:nbStates), "total")
    lty <- c(rep(1, nbStates), 2)
    lwd <- c(rep(1, nbStates), 2)
    lineCol <- c(col, "black")

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
        lines(distData$step$step, distData$step[,i+1], col = lineCol[i],
              lty = lty[i], lwd = lwd[i])
    }
    legend("top", legText, lwd = lwd, col = lineCol, lty = lty, bty = "n")

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
        lines(distData$angle$angle, distData$angle[,i+1], col = lineCol[i],
              lty = lty[i], lwd = lwd[i])
    }
    legend("top", legText, lwd = lwd, col = lineCol, lty = lty, bty = "n")

    ##################################################
    ## Plot the t.p. as functions of the covariates ##
    ##################################################
    beta <- m$mle$beta
    if(nbStates > 1 & nrow(beta) > 1) {
        trProbs <- getPlotData(m, type = "tpm", format = "wide")

        # loop over covariates
        par(mfrow = c(nbStates, nbStates))
        par(mar = c(5, 4, 4, 2) - c(0, 0, 1.5, 1)) # bottom, left, top, right
        for(cov in 1:ncol(m$rawCovs)) {
            trProbsCov <- trProbs[[cov]]
            # loop over entries of the transition probability matrix
            for(i in 1:nbStates) {
                for(j in 1:nbStates) {
                    trName <- paste0("S", i, "toS", j)
                    plot(trProbsCov[,1], trProbsCov[,trName], type = "l",
                         ylim = c(0, 1), xlab = names(trProbs)[cov],
                         ylab = paste(i, "->", j))

                    # derive confidence intervals using the delta method
                    if(plotCI) {
                        options(warn = -1) # to muffle "zero-length arrow..." warning
                        # plot the confidence intervals
                        arrows(trProbsCov[,1], trProbsCov[,paste0(trName, ".lci")],
                               trProbsCov[,1], trProbsCov[,paste0(trName, ".uci")],
                               length = 0.025, angle = 90, code = 3,
                               col = gray(0.5), lwd = 0.7)
                        options(warn = 1)
                    }
                }
            }

            mtext("Transition probabilities", side = 3, outer = TRUE, padj = 2)
        }
    }

    #################################
    ## Plot maps colored by states ##
    #################################
    # Prepare the data
    nbAnimals <- length(unique(m$data$ID))
    if(is.character(animals)) {
        if(any(!animals %in% unique(m$data$ID))) {
            stop("Check 'animals' argument, ID not found")
        }
        animalsInd <- which(unique(m$data$ID) %in% animals)
    } else if(is.numeric(animals)) {
        if(min(animals) < 1 | max(animals) > nbAnimals) {
            stop("Check 'animals' argument, index out of bounds")
        }
        animalsInd <- animals
    } else {
        animalsInd <- 1:nbAnimals
    }
    nbAnimals <- length(animalsInd)
    ID <- unique(m$data$ID)[animalsInd]

    if(nbStates>1 & plotTracks) { # no need to plot the map if only one state
        par(mfrow = c(1, 1))
        par(mar = c(5, 4, 4, 2) - c(0, 0, 2, 1)) # bottom, left, top, right

        for(zoo in 1:nbAnimals) {
            # data for this animal
            ind <- which(m$data$ID == ID[zoo])
            s <- states[ind]
            x <- m$data$x[ind]
            y <- m$data$y[ind]

            # slightly different for 2D and 1D data
            if(!all(y == 0)) {
                plot(x, y, pch = 16, col = col[s], cex = 0.5, asp = 1,
                     xlab = "x", ylab = "y")
                segments(x0 = x[-length(x)], y0 = y[-length(y)],
                         x1 = x[-1], y1 = y[-1],
                         col = col[s[-length(s)]], lwd = 1.3)
            } else { # if 1D data
                plot(x, xlab = "time", ylab = "x", pch = 16,
                     cex = 0.5, col = col[s])
                segments(x0 = 1:(length(x) - 1), y0 = x[-length(x)],
                         x1 = 2:length(x), y1 = x[-1],
                         col = col[s[-length(x)]], lwd = 1.3)
            }

            mtext(paste("Animal ID:", ID[zoo]), side = 3, outer = TRUE, padj = 2)
        }
    }

    # set the graphical parameters back to default
    par(mfrow = c(1, 1))
    par(mar = c(5, 4, 4, 2) + 0.1) # bottom, left, top, right
    par(ask = FALSE)
}
