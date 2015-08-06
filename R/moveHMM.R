
#' Constructor of moveHMM objects
#'
#' @param mod A list of attributes of the fitted model. Mandatory slots : data, states, mle,
#' stepDist, angleDist, mod (result of call to optimizer),
#' Optional slots : angleMean (if not estimated)
#'
#' @return An object moveHMM.

moveHMM <- function(m)
{
  if(is.null(m$data) | is.null(m$mle) | is.null(m$stepDist) |
       is.null(m$angleDist) | is.null(m$mod))
    stop("Can't construct moveHMM object : fields are missing")

  obj <- m

  class(obj) <- append(class(obj),"moveHMM")
  return(obj)
}
