
#' Pseudo-residuals
#'
#' The pseudo-residuals of a moveHMM model, as described in Zucchini and McDonad (2009).
#'
#' @param m A \code{moveHMM} object.
#'
#' @return A list of:
#' \item{stepRes}{The pseudo-residuals for the step lengths}
#' \item{angleRes}{The pseudo-residuals for the turning angles}
#'
#' @details If some turning angles in the data are equal to pi, the corresponding pseudo-residuals
#' will not be included. Indeed, given that the turning angles are defined on (-pi,pi], an angle of pi
#' results in a pseudo-residual of +Inf (check Section 6.2 of reference for more information on the
#' computation of pseudo-residuals).
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#' res <- pseudoRes(m)
#' qqnorm(res$stepRes)
#' qqnorm(res$angleRes)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).
#'
#' @export
#' @importFrom stats integrate qnorm

pseudoRes <- function(m)
{
    if(!is.moveHMM(m))
        stop("'m' must be a moveHMM object (as output by fitHMM)")

    stepFun <- paste("p",m$conditions$stepDist,sep="")

    angleDist <- m$conditions$angleDist
    if(angleDist!="none") {
        angleFun <- paste("d",angleDist,sep="") # integrated below

        if(length(which(m$data$angle==pi))>0)
            message("Note: Some angles are equal to pi, and the corresponding pseudo-residuals are not included")
    }

    data <- m$data
    nbObs <- nrow(data)
    nbStates <- ncol(m$mle$stepPar)

    # forward log-probabilities
    la <- logAlpha(m)

    stepRes <- rep(NA,nbObs)
    pStepMat <- matrix(NA,nbObs,nbStates)

    if(angleDist!="none") {
        angleRes <- rep(NA,nbObs)
        pAngleMat <- matrix(NA,nbObs,nbStates)
    }
    else {
        angleRes <- NULL
        pAngleMat <- NULL
    }

    for(state in 1:nbStates) {
        # define lists of parameters
        stepArgs <- list(data$step[1])
        if(!m$conditions$zeroInflation) {
            for(k in 1:nrow(m$mle$stepPar))
                stepArgs[[k+1]] <- m$mle$stepPar[k,state]

            zeromass <- 0
        }
        else {
            for(k in 1:(nrow(m$mle$stepPar)-1))
                stepArgs[[k+1]] <- m$mle$stepPar[k,state]

            zeromass <- m$mle$stepPar[nrow(m$mle$stepPar),state]
        }

        if(m$conditions$stepDist=="gamma") {
            shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
            scale <- stepArgs[[3]]^2/stepArgs[[2]]
            stepArgs[[2]] <- shape
            stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
        }

        if(angleDist!="none") {
            angleArgs <- list(angleFun,-pi,data$angle[1]) # to pass to function "integrate" below
            for(k in 1:nrow(m$mle$anglePar))
                angleArgs[[k+3]] <- m$mle$anglePar[k,state]
        }

        for(i in 1:nbObs) {
            if(!is.na(data$step[i])) {
                stepArgs[[1]] <- data$step[i]
                pStepMat[i,state] <- zeromass+(1-zeromass)*do.call(stepFun,stepArgs)
            }

            if(angleDist!="none") {
                if(!is.na(data$angle[i])) {
                    # angle==pi => residual=Inf
                    if(data$angle[i]!=pi) {
                        angleArgs[[3]] <- data$angle[i]
                        pAngleMat[i,state] <- do.call(integrate,angleArgs)$value
                    }
                }
            }
        }
    }

    # define covariates
    covs <- model.matrix(m$conditions$formula, data = data)

    if(nbStates>1)
        trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))
    else
        trMat <- array(1,dim=c(1,1,nbObs))

    # aInd = list of indices of first observation for each animal
    aInd <- c(1, which(data$ID[-1] != data$ID[-nbObs]) + 1)

    for(i in 1:nbObs) {
        if(!is.na(data$step[i])){
            if(any(i==aInd)){

                stepRes[i] <- qnorm(t(m$mle$delta)%*%pStepMat[i,])

                if(angleDist!="none") {
                    if(!is.na(data$angle[i]))
                        angleRes[i] <- qnorm(t(m$mle$delta)%*%pAngleMat[i,])
                }
            } else {
                gamma <- trMat[,,i]
                c <- max(la[i-1,]) # cancels below ; prevents numerical errors
                a <- exp(la[i-1,]-c)

                stepRes[i] <-qnorm(t(a)%*%(gamma/sum(a))%*%pStepMat[i,])

                if(angleDist!="none") {
                    if(!is.na(data$angle[i]))
                        angleRes[i] <- qnorm(t(a)%*%(gamma/sum(a))%*%pAngleMat[i,])
                }
            }
        }
    }

    return(list(stepRes=stepRes,angleRes=angleRes))
}
