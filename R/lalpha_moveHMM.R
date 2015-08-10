
#' Forward probabilities
#'
#' @param m A moveHMM object.
#'
#' @return The matrix of forward probabilities.

lalpha <- function(m) UseMethod("lalpha") # define generic method lalpha

lalpha.moveHMM <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  lalpha <- matrix(NA,nbObs,nbStates)

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- data[,covsCol]

  allProbs <- allProbs(data,nbStates,m$stepDist,m$angleDist,m$mle$stepPar,m$mle$anglePar,m$zeroInflation)
  trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))

  lscale <- 0
  foo <- m$mle$delta*allProbs[1,]
  lalpha[1,] <- log(foo)+lscale

  for(i in 2:nbObs) {
    gamma <- trMat[,,i]
    foo <- foo%*%gamma*allProbs[i,]
    lscale <- lscale+log(sum(foo))
    foo <- foo/sum(foo)
    lalpha[i,] <- log(foo)+lscale
  }

  return(lalpha)
}
