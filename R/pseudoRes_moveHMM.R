
#' Generic pseudoRes method
#' @param m Fitted model
pseudoRes <- function(m) UseMethod("pseudoRes") # define generic method pseudoRes

#' Pseudo-residuals
#' @method pseudoRes moveHMM
#'
#' @param m A \code{moveHMM} object.
#'
#' @return A list of :
#' \item{stepRes}{The pseudo-residuals for the step lengths}
#' \item{angleRes}{The pseudo-residuals for the turning angles}
#'
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#' res <- pseudoRes(m)
#' qqnorm(res$stepRes)
#' qqnorm(res$angleRes)
#'
#' @references
#' Zucchini, W. and MacDonald, I.L. 2009.
#' Hidden Markov Models for Time Series: An Introduction Using R.
#' Chapman & Hall (London).

pseudoRes.moveHMM <- function(m)
{
  stepFun <- paste("p",m$stepDist,sep="")

  angleDist <- m$angleDist
  if(angleDist!="none") {
    angleFun <- paste("d",angleDist,sep="") # integrated below

    if(length(which(data$angle==pi))>0)
      warning("Some angle are equal to pi, and the corresponding pseudo-residuals are not included")
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

    if(m$stepDist=="gamma") {
      shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
      scale <- stepArgs[[3]]^2/stepArgs[[2]]
      stepArgs[[2]] <- shape
      stepArgs[[3]] <- 1/scale # dgamma expects rate=1/scale
    }

    if(angleDist!="none") {
      angleArgs <- list(angleFun,-pi,data$angle[1]) # to pass to function "integrate" below
      for(k in 1:nrow(m$mle$anglePar))
        angleArgs[[k+3]] <- m$mle$anglePar[k,state]

      for(i in 1:nbObs) {
        if(!is.na(data$step[i])) {
          stepArgs[[1]] <- data$step[i]
          pStepMat[i,state] <- zeromass+(1-zeromass)*do.call(stepFun,stepArgs)
        }

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

  if(!is.na(data$step[1]))
    stepRes[1] <- qnorm(t(m$mle$delta)%*%pStepMat[1,])

  if(angleDist!="none") {
    if(!is.na(data$angle[1]))
      angleRes[1] <- qnorm(t(m$mle$delta)%*%pAngleMat[1,])
  }

  # define covariates
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- data[,covsCol]

  trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))

  for(i in 2:nbObs) {
    gamma <- trMat[,,i]
    c <- max(la[i-1,]) # cancels below ; prevents numerical errors
    a <- exp(la[i-1,]-c)

    if(!is.na(data$step[i]))
      stepRes[i] <-qnorm(t(a)%*%(gamma/sum(a))%*%pStepMat[i,])

    if(angleDist!="none") {
      if(!is.na(data$angle[i]))
        angleRes[i] <- qnorm(t(a)%*%(gamma/sum(a))%*%pAngleMat[i,])
    }
  }

  return(list(stepRes=stepRes,angleRes=angleRes))
}
