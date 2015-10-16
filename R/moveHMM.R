
#' Constructor of \code{moveHMM} objects
#'
#' @param m A list of attributes of the fitted model: \code{mle} (the maximum likelihood estimates of
#' the parameters of the model), \code{data} (the movement data), \code{stepDist} (the step length
#' distribution name), \code{angleDist} (the turning angle distribution name), \code{mod} (the object
#' returned by the numerical optimizer \code{nlm}), \code{conditions} (a few conditions used to fit
#' the model: \code{zeroInflation}, \code{estAngleMean}, \code{stationary}, and \code{formula}),
#' \code{rawCovs} (optional -- only if there are covariates in the data).
#'
#' @return An object \code{moveHMM}.
#'
#' @export

moveHMM <- function(m)
{
  if(is.null(m$data) | is.null(m$mle) | is.null(m$stepDist) |
     is.null(m$angleDist) | is.null(m$mod) | is.null(m$conditions))
    stop("Can't construct moveHMM object: fields are missing")

  obj <- m

  class(obj) <- append("moveHMM",class(obj))
  return(obj)
}
