\name{trMatrix_rcpp}
\alias{trMatrix_rcpp}
\title{
    Transition probability matrix
}
\description{
    Computation of the transition probability matrix, as a function of the covariate values and the regression parameters beta. Written in C++. Used in function viterbi.
}
\arguments{
    \item{nbStates}{Number of states.}
    \item{beta}{Matrix of regression parameters.}
    \item{covs}{Matrix of covariate values}
}
\value{
    Three-dimensional array trMat, such that trMat[,,t] is the transition probability matrix at time t.
}
