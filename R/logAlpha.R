
#' Forward log-probabilities
#'
#' Used in \code{\link{stateProbs}} and \code{\link{pseudoRes}}.
#'
#' @param m A \code{moveHMM} object.
#'
#' @return The matrix of forward log-probabilities.
#'
#' @examples
#' \dontrun{
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' la <- logAlpha(m)
#' }

logAlpha <- function(m)
{
  data <- m$data
  nbStates <- ncol(m$mle$stepPar)
  nbObs <- nrow(data)
  lalpha <- matrix(NA,nbObs,nbStates)

  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle")
  covs <- data[,covsCol]

  allProbs <- allProbs(data,nbStates,m$conditions$stepDist,m$conditions$angleDist,m$mle$stepPar,
                       m$mle$anglePar,m$conditions$zeroInflation)

  if(nbStates>1)
    trMat <- trMatrix_rcpp(nbStates,m$mle$beta,as.matrix(covs))
  else
    trMat <- array(1,dim=c(1,1,nbObs))

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
