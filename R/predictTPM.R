
#' Predict transition probabilities for new covariate values
#'
#' @param m Fitted moveHMM object, as returned by \code{\link{fitHMM}}
#' @param newData Data frame with columns for the covariates
#' @param beta Optional matrix of regression coefficients for the transition
#' probability model. By default, uses estimates in \code{m}.
#' @param returnCI Logical indicating whether confidence intervals should
#' be returned. Default: FALSE.
#' @param alpha Confidence level if returnCI = TRUE. Default: 0.95, i.e.,
#' 95\% confidence intervals.
#'
#' @return List with elements 'mle', 'lci', and 'uci' (the last two only if
#' returnCI = TRUE). Each element is an array, where each layer is a
#' transition probability matrix corresponding to a row of newData.
#'
#' @export
predictTPM <- function(m, newData, beta = m$mle$beta,
                       returnCI = FALSE, alpha = 0.95) {
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    if(length(m$mod)<=1)
        stop("The given model hasn't been fitted.")

    if(alpha<0 | alpha>1)
        stop("alpha needs to be between 0 and 1.")

    nbStates <- ncol(m$mle$stepPar)

    if(nbStates>1) {
        # covariance matrix of estimates
        if(!is.null(m$mod$hessian)) {
            Sigma <- ginv(m$mod$hessian)
        } else {
            returnCI <- FALSE
        }

        # indices corresponding to regression coefficients in m$mod$estimate
        i1 <- length(m$mle$stepPar) + length(m$mle$anglePar) -
            (!m$conditions$estAngleMean)*nbStates + 1
        i2 <- i1 + length(beta) - 1
        gamInd <- i1:i2

        quantSup <- qnorm(1 - (1 - alpha) / 2)

        # Make design matrix for new data
        desMat <- model.matrix(m$conditions$formula,data = newData)

        # Get transition probabilities
        trMat <- trMatrix_rcpp(nbStates, beta, desMat)
        rownames(trMat) <- paste0("fromState", 1:nbStates)
        colnames(trMat) <- paste0("toState", 1:nbStates)
        out <- list(mle = trMat)

        # loop over entries of the transition probability matrix
        if(returnCI) {
            # for differentiation to obtain confidence intervals (delta method)
            get_gamma <- function(beta, covs, nbStates, i, j){
                gamma <- trMatrix_rcpp(nbStates, beta, covs)[,,1]
                return(gamma[i, j])
            }

            # Copy structure of trMat for upper/lower CI bounds
            lci <- trMat
            uci <- trMat

            # Loop over transition probabilities
            for(i in 1:nbStates) {
                for(j in 1:nbStates) {
                    # derive confidence intervals using the delta method
                    dN <- t(apply(desMat, 1, function(x)
                        grad(get_gamma, beta, covs = matrix(x, nrow = 1),
                                       nbStates = nbStates, i = i,j = j)))

                    se <- t(apply(dN, 1, function(x)
                        suppressWarnings(sqrt(x %*% Sigma[gamInd, gamInd] %*% x))))

                    # transform estimates and standard errors to R, to derive CI on working scale,
                    # then back-transform to [0,1]
                    width <- quantSup * se / (trMat[i,j,] - trMat[i,j,]^2)
                    lci[i,j,] <- plogis(qlogis(trMat[i,j,]) - width)
                    uci[i,j,] <- plogis(qlogis(trMat[i,j,]) + width)
                }
            }

            out$lci <- lci
            out$uci <- uci
        }
    } else {
        stop("Only 1 state -- no transition probabilities to predict.")
    }

    return(out)
}
