
#' Constructor of \code{moveHMM} objects
#'
#' This constructor is not intended to be used, except inside the function
#' \code{\link{fitHMM}}. Refer to the documentation for that function.
#'
#' @param m A list of attributes of the fitted model: \code{mle} (maximum
#' likelihood estimates of model parameters), \code{data} (movement data),
#' \code{mod} (output of numerical optimizer \code{nlm}), \code{conditions}
#' (few conditions used to fit the model: \code{stepDist}, \code{angleDist},
#' \code{zeroInflation}, \code{estAngleMean}, \code{stationary}, and
#' \code{formula}), \code{rawCovs} (optional -- only if there are covariates
#' in the data).
#'
#' @return An object \code{moveHMM}.

moveHMM <- function(m)
{
    if(is.null(m$data) | is.null(m$mle) | is.null(m$mod) | is.null(m$conditions))
        stop("Can't construct moveHMM object: fields are missing")

    obj <- m

    class(obj) <- append("moveHMM",class(obj))
    return(obj)
}

#' Is moveHMM
#'
#' Check that an object is of class \code{\link{moveHMM}}. Used in \code{\link{CI}},
#' \code{\link{plotPR}}, \code{\link{plotStates}}, \code{\link{pseudoRes}}, \code{\link{stateProbs}},
#' and \code{\link{viterbi}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{moveHMM}}, \code{FALSE} otherwise.

is.moveHMM <- function(x)
    inherits(x,"moveHMM")
