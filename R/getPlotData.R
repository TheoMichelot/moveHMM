
#' Data to produce plots of fitted model
#'
#' @param m Fitted HMM object, as output by fitHMM.
#' @param type Type of plot, one of: "dist", "tpm", "stat"
#' @param format Format of data, either "wide" (for base graphics)
#' or "long" (for ggplot)
#' @param alpha Level of confidence intervals. Default: 0.95, i.e.,
#' 95\% confidence intervals
#'
#' @details
#' \itemize{
#' \item{If type = "dist", the function evaluates each state-dependent
#' distribution over the range of observed variable (step length or
#' turning angle), and weighs them by the proportion of time spent
#' in each state (obtained from Viterbi state sequence).}
#' \item{If type = "tpm", the function returns transition probabilities
#' estimated over a range of covariate values. Other covariates are
#' fixed to their mean values.}
#' }
#'
#' @return Data frame (or list of data frames) containing data in a format
#' that can easily be plotted. If type = "dist", the output is a list with
#' two elements, "step" and "angle". If type = "tpm" or "stat", the output
#' is a list with one element for each covariate. See details for more
#' extensive description of output.
#'
#' @export
getPlotData <- function(m, type, format = "wide", alpha = 0.95) {
    nbStates <- ncol(m$mle$stepPar)
    out <- list()

    if(type == "dist") {
        ###################################
        ## State-dependent distributions ##
        ###################################
        # proportion of time steps in each state
        w <- 1
        if(nbStates > 1) {
            states <- viterbi(m)
            w <- sapply(1:3, function(s) length(which(states == s))/length(states))
        }

        ##################
        ## Step lengths ##
        ##################
        # unpack step length parameters
        if(m$conditions$zeroInflation) {
            zeromass <- m$mle$stepPar[nrow(m$mle$stepPar),]
            stepPar <- as.matrix(m$mle$stepPar[-nrow(m$mle$stepPar),])
        } else {
            stepPar <- m$mle$stepPar
        }

        stepFun <- paste0("d", m$conditions$stepDist)
        stepGrid <- seq(0, max(m$data$step, na.rm = TRUE), length = 1000)
        stepDensities <- data.frame(step = stepGrid)
        for(state in 1:nbStates) {
            stepArgs <- list(stepGrid)
            for(j in 1:nrow(stepPar)) {
                stepArgs[[j+1]] <- stepPar[j,state]
            }

            # conversion between mean/sd and shape/rate if necessary
            if(m$conditions$stepDist == "gamma") {
                shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
                rate <- stepArgs[[2]]/stepArgs[[3]]^2
                stepArgs[[2]] <- shape
                stepArgs[[3]] <- rate
            }
            # (weighted by the proportion of each state in the Viterbi states sequence)
            stepDensities[[paste0("state", state)]] <- w[state] * do.call(stepFun,stepArgs)
            if(m$conditions$zeroInflation) {
                stepDensities[[paste0("state", state)]] <-
                    (1-zeromass[state]) * stepDensities[[paste0("state", state)]]
            }
        }
        stepDensities$total <- rowSums(stepDensities[,-1])

        # wide or long format
        out$step <- stepDensities
        if(format == "long") {
            out$step <- data.frame(
                step = out$step[,1],
                density = unlist(out$step[,-1], use.names = FALSE),
                state = rep(c(paste0("state ", 1:nbStates), "total"), each = nrow(out$step)))
        }

        ####################
        ## Turning angles ##
        ####################
        angleDensities <- NULL
        if(m$conditions$angleDist != "none") {
            angleFun <- paste("d", m$conditions$angleDist, sep = "")
            angleGrid <- seq(-pi, pi, length = 1000)
            angleDensities <- data.frame(angle = angleGrid)

            for(state in 1:nbStates) {
                angleArgs <- list(angleGrid)
                for(j in 1:nrow(m$mle$anglePar)) {
                    angleArgs[[j+1]] <- m$mle$anglePar[j,state]
                }

                # (weighted by the proportion of each state in the Viterbi states sequence)
                angleDensities[[paste0("state", state)]] <- w[state] * do.call(angleFun, angleArgs)
            }
            angleDensities$total <- rowSums(angleDensities[,-1])
        }

        # wide or long format
        out$angle <- angleDensities
        if(format == "long") {
            out$angle <- data.frame(
                angle = out$angle[,1],
                density = unlist(out$angle[,-1], use.names = FALSE),
                state = rep(c(paste0("state ", 1:nbStates), "total"), each = nrow(out$angle)))
        }
    } else if(type == "tpm" | type == "stat") {
        beta <- m$mle$beta
        if(nrow(beta) == 1) {
            stop("No covariate effects on tpm.")
        }
        rawCovs <- m$rawCovs
        gridLength <- 100

        # loop over covariates
        for(cov in 1:ncol(m$rawCovs)) {
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

            if(type == "tpm") {
                ##############################
                ## Transition probabilities ##
                ##############################
                # Transition probabilities on covariate grid
                trMat <- predictTPM(m = m, newData = tempCovs,
                                    returnCI = TRUE, alpha = 0.95)

                # all transition probabilities in the right order for formatting below
                trMatVec <- unlist(lapply(trMat, function(a) aperm(a, c(3, 2, 1))))

                if(format == "wide") {
                    # wide format: one column for covariate, one column for each t.p.,
                    # one column for each lower bound, one column for each upper bound
                    plotData <- matrix(c(tempCovs[,cov], trMatVec),
                                       nrow = nrow(tempCovs))
                    colnames(plotData) <- c(
                        colnames(rawCovs)[cov],
                        paste0("S", rep(rep(1:nbStates, each = nbStates), 3),
                               "toS", rep(rep(1:nbStates, nbStates), 3),
                               rep(c("", ".lci", ".uci"), each = nbStates^2)))
                    plotData <- as.data.frame(plotData)
                } else if(format == "long") {
                    # long format: one column for covariate, one column for mle,
                    # one column for lower bound, one column for upper bound, and
                    # one column for t.p. entry name
                    plotData <- as.data.frame(matrix(
                        c(rep(tempCovs[,cov], nbStates^2), trMatVec), ncol = 4))
                    colnames(plotData) <- c(colnames(rawCovs)[cov], "mle", "lci", "uci")
                    plotData$prob <- rep(paste0("Pr(", rep(1:nbStates, each = nbStates),
                                                " -> ", rep(1:nbStates, nbStates), ")"),
                                         each = nrow(tempCovs))
                }
            } else if(type == "stat") {
                ####################################
                ## Stationary state probabilities ##
                ####################################
                # State probabilities on covariate grid
                stat <- predictStationary(m = m, newData = tempCovs,
                                          returnCI = TRUE, alpha = 0.95)
                statVec <- unlist(stat, use.names = FALSE)

                if(format == "wide") {
                    # wide format: one column for covariate, one column for each prob,
                    # one column for each lower bound, one column for each upper bound
                    plotData <- matrix(c(tempCovs[,cov], statVec),
                                       nrow = nrow(tempCovs))
                    colnames(plotData) <- c(
                        colnames(rawCovs)[cov],
                        paste0("S", rep(1:nbStates, 3),
                               rep(c("", ".lci", ".uci"), each = nbStates)))
                    plotData <- as.data.frame(plotData)
                } else if(format == "long") {
                    # long format: one column for covariate, one column for mle,
                    # one column for lower bound, one column for upper bound, and
                    # one column for state
                    plotData <- as.data.frame(matrix(
                        c(rep(tempCovs[,cov], nbStates), statVec), ncol = 4))
                    colnames(plotData) <- c(colnames(rawCovs)[cov], "mle", "lci", "uci")
                    plotData$state <- rep(1:nbStates, each = nrow(tempCovs))
                }
            }

            out[[colnames(rawCovs)[cov]]] <- plotData
        }
    }

    return(out)
}
