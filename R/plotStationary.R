
#' Plot stationary state probabilities
#'
#' @param m An object \code{moveHMM}
#' @param col Vector or colors for the states (one color per state).
#' @param plotCI Logical. Should 95\% confidence intervals be plotted? (Default: FALSE)
#' @param alpha Significance level of the confidence intervals if plotCI=TRUE.
#' Default: 0.95 (i.e. 95\% CIs).
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' plotStationary(m)
#'
#' @export
plotStationary <- function(m, col=NULL, plotCI=FALSE, alpha=0.95)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    data <- m$data
    nbStates <- ncol(m$mle$stepPar)
    beta <- m$mle$beta

    if(nrow(beta)==1)
        stop("No covariate effect to plot (nrow(beta)==1).")

    # prepare colors for the states
    if(is.null(col) | (!is.null(col) & length(col) != nbStates)) {
        col <- getPalette(nbStates = nbStates)
    }

    rawCovs <- m$rawCovs
    gridLength <- 100

    # loop over covariates
    for(cov in 1:ncol(rawCovs)) {
        inf <- min(rawCovs[,cov], na.rm = TRUE)
        sup <- max(rawCovs[,cov], na.rm = TRUE)

        # mean values of each covariate
        meanCovs <- colMeans(rawCovs)

        # set all covariates to their mean, except for "cov"
        # (which takes a grid of values from inf to sup)
        tempCovs <- as.data.frame(matrix(rep(meanCovs, each = gridLength),
                                         nrow = gridLength))
        tempCovs[,cov] <- seq(inf, sup, length = gridLength)
        colnames(tempCovs) <- colnames(rawCovs)

        stat <- predictStationary(m = m, newData = tempCovs,
                                   returnCI = TRUE, alpha = 0.95)

        # plot lines
        plot(NA, xlim = c(inf, sup), ylim = c(0, 1), xlab = names(rawCovs)[cov],
             ylab = "Stationary state probabilities")
        for(state in 1:nbStates) {
            lines(tempCovs[,cov], stat$mle[,state], col = col[state], lwd = 2)
        }
        legend("topleft", legend = paste("state", 1:nbStates), col = col,
               lty = 1, bty = "n")

        # add confidence intervals
        if(plotCI) {
            for(state in 1:nbStates) {
                options(warn = -1) # to muffle "zero-length arrow..." warning
                # plot the confidence intervals
                arrows(tempCovs[,cov], stat$lci[,state],
                       tempCovs[,cov], stat$uci[,state],
                       length = 0.025, angle = 90, code = 3,
                       col = col[state], lwd = 0.7)
                options(warn = 1)
            }
        }
    }
}
