
#' Plot stationary state probabilities
#'
#' @param m An object \code{moveHMM}
#' @param col Vector or colors for the states (one color per state).
#' @param plotCI Logical. Should 95\% confidence intervals be plotted? (Default: FALSE)
#'
#' @export
plotStationary <- function(m, col=NULL, plotCI=FALSE)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    data <- m$data
    nbStates <- ncol(m$mle$stepPar)
    beta <- m$mle$beta

    if(nrow(beta)==1)
        stop("No covariate effect to plot (nrow(beta)==1).")

    # prepare colors for the states
    if(!is.null(col) & length(col)!=nbStates) {
        warning("Length of 'col' should be equal to number of states - argument ignored")
        col <- NULL
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

    # for differentiation in delta method below
    get_stat <- function(beta,covs,nbStates,i) {
        gamma <- momentuHMM:::trMatrix_rcpp(nbStates,beta,covs)[,,1]
        solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i]
    }

    rawCovs <- m$rawCovs
    gridLength <- 100

    # loop over covariates
    for(cov in 1:ncol(rawCovs)) {
        inf <- min(rawCovs[,cov],na.rm=TRUE)
        sup <- max(rawCovs[,cov],na.rm=TRUE)

        # mean values of each covariate
        meanCovs <- colMeans(rawCovs)

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
            probs <- stationary(m, covMat=desMat)

            plot(tempCovs[,cov], probs[,1], type="l", ylim=c(0,1), col=col[1],
                 xlab=names(rawCovs)[cov], ylab="Stationary state probabilities")
            for(state in 2:nbStates)
                points(tempCovs[,cov], probs[,state], type="l", col=col[state])
            legend("topleft", legend=paste("State",1:nbStates), col=col, lty=1, bty="n")

            if(plotCI) {
                # covariance matrix of estimates
                Sigma <- ginv(m$mod$hessian)

                # indices corresponding to regression coefficients in m$mod$estimate
                i1 <- length(m$mle$stepPar) + length(m$mle$anglePar) + 1
                i2 <- i1 + length(m$mle$beta) - 1
                gamInd <- i1:i2

                lci <- matrix(NA,gridLength,nbStates)
                uci <- matrix(NA,gridLength,nbStates)

                for(state in 1:nbStates) {
                    dN <- t(apply(desMat, 1, function(x)
                        numDeriv::grad(get_stat,beta,covs=matrix(x,nrow=1),nbStates=nbStates,i=state)))

                    se <- t(apply(dN, 1, function(x)
                        suppressWarnings(sqrt(x%*%Sigma[gamInd,gamInd]%*%x))))

                    lci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) -
                                                    1.96 * (1/(probs[,state]-probs[,state]^2)) * se)))
                    uci[,state] <- 1/(1 + exp(-(log(probs[,state]/(1-probs[,state])) +
                                                    1.96 * (1/(probs[,state]-probs[,state]^2)) * se)))

                    # plot the confidence intervals
                    arrows(tempCovs[,cov], lci[,state], tempCovs[,cov], uci[,state], length=0.025,
                           angle=90, code=3, col=col[state], lwd=0.7)
                }
            }
        }
    }
}
