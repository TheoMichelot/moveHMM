
#' Generate the transition probability matrix
#'
#' Used in the computation of the log-likelihood.
#'
#' @param nbStates Number of states of the HMM.
#' @param beta Matrix of regression parameters. Number of rows : number of covariates + 1 ;
#' number of columns : number of non-diagonal elements in the t.p.m. (i.e. nbStates*(nbStates-1)).
#' @param covs Design matrix of covariates.
#'
#' @return A three-dimensional array gamma, such that gamma[,,t] is the transition probability
#' matrix corresponding to observation t.

trMatrix <- function(nbStates,beta,covs)
{
  nbObs <- nrow(covs)
  gamma <- array(0,c(nbStates,nbStates,nbObs))
  for(i in 1:nbStates) gamma[i,i,] <- 1 # diagonals of one on each layer

  g <- exp(as.matrix(covs)%*%beta)
  gamma[!gamma] <- t(g) # transpose because R picks elements column-wise

  gamma <- aperm(gamma,c(2,1,3)) # transpose because R fills elements colum-wise
  s <- array(NA,c(nbStates,nbStates,nbObs))
  for(i in 1:nbStates) s[,i,] <- apply(gamma,c(1,3),sum) # row sums
  gamma <- gamma/s

  return(gamma)
}
