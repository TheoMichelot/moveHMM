
#' Constructor of moveHMM objects
#'
#' @param m A list of attributes of the fitted model. Mandatory slots : data, states, mle,
#' stepDist, angleDist, mod (result of call to optimizer), conditions (zeroInflation,
#' estAngleMean, and stationary).
#'
#' @return An object moveHMM.

moveHMM <- function(m)
{
  if(is.null(m$data) | is.null(m$mle) | is.null(m$stepDist) |
     is.null(m$angleDist) | is.null(m$mod) | is.null(m$states) |
     is.null(m$conditions))
    stop("Can't construct moveHMM object : fields are missing")

  obj <- m

  class(obj) <- append("moveHMM",class(obj))
  return(obj)
}
