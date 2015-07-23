\name{nLogLike_rcpp}
\alias{nLogLike_rcpp}
\title{
    Forward algorithm
}
\description{
    Forward algorithm computation. This is the main loop of the likelihood computation, written in C++.
}
\arguments{
    \item{trVec}{Array of the transition probabilities. Passed as an array to C++, and then built back into a 3-dimensional array.}
    \item{delta}{Stationary distribution.}
    \item{allProbs}{Matrix of joint probabilities.}
}
\value{
    Log-likelihood value.
}
