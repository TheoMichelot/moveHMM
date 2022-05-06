
#' Predict stationary state probabilities
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
#' returnCI = TRUE). Each element is a matrix of stationary state probabilities
#' with one row for each row of newData and one column for each state.
#'
#' @export
predictStationary <- function(m, newData, beta = m$mle$beta,
                              returnCI = FALSE, alpha = 0.95) {
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    data <- m$data
    nbStates <- ncol(m$mle$stepPar)

    # for differentiation in delta method below
    get_stat <- function(beta,covs,nbStates,i) {
        gamma <- trMatrix_rcpp(nbStates,beta,covs)[,,1]
        solve(t(diag(nbStates)-gamma+1),rep(1,nbStates))[i]
    }

    # for confidence intervals
    quantSup <- qnorm(1-(1-alpha)/2)

    # stationary probs
    desMat <- model.matrix(m$conditions$formula, data = newData)
    probs <- stationary(m, covs = desMat, beta = beta)
    colnames(probs) <- paste0("state", 1:nbStates)
    out <- list(mle = probs)

    if(returnCI) {
        # copy structure of probs matrix
        lci <- probs
        uci <- probs

        # covariance matrix of estimates
        Sigma <- ginv(m$mod$hessian)

        # indices corresponding to regression coefficients in m$mod$estimate
        i1 <- length(m$mle$stepPar) + length(m$mle$anglePar) -
            (!m$conditions$estAngleMean)*nbStates + 1
        i2 <- i1 + length(m$mle$beta) - 1
        gamInd <- i1:i2

        for(state in 1:nbStates) {
            dN <- t(apply(desMat, 1, function(x)
                grad(get_stat, beta, covs = matrix(x, nrow = 1),
                     nbStates = nbStates, i = state)))

            se <- t(apply(dN, 1, function(x)
                suppressWarnings(sqrt(x %*% Sigma[gamInd, gamInd] %*% x))))

            # transform estimates and standard errors to R, to derive CI on working scale,
            # then back-transform to [0,1]
            width <- quantSup*se/(probs[,state]-probs[,state]^2)
            lci[,state] <- plogis(qlogis(probs[,state]) - width)
            uci[,state] <- plogis(qlogis(probs[,state]) + width)
        }

        out$lci <- lci
        out$uci <- uci
    }

    return(out)
}
