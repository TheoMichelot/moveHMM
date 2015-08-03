\name{nLogLike_rcpp}
\alias{nLogLike_rcpp}
\title{
    Negative log-likelihood function
}
\description{
    Computation of the negative log-likelihood. This function computes the transition probability matrix and the matrix of joint probablities, and then runs the forward algorithm to compute the log-likelihood. Written in C++.
}
\arguments{
    \item{nbStates}{Number of states}
    \item{beta}{Matrix of regression coefficients for the transition probability matrix.}
    \item{covs}{Covariates}
    \item{data}{The observations, in a dataframe.}
    \item{stepDist}{The name of the step length distribution.}
    \item{angleDist}{The name of the turning angle distribution (defaults to "NULL").}
    \item{stepPar}{State-dependent parameters of the step length distribution.}
    \item{anglePar}{State-dependent parameters of the turning angle distribution (defaults to NULL).}
    \item{delta}{Stationary distribution.}
    \item{aInd}{Vector of indices of the rows at which the data switches to another animal.}
    \item{zeroInflation}{true if zero inflation is included in the step length distribution, false otherwise (default).}
}
\value{
    Negative log-likelihood.
}
