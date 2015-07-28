
#' Constructor of moveHMM objects
#'
#' @param mod A list of attributes of the fitted model. Mandatory slots : data, states, mle,
#' stepDist, angleDist,
#' Optional slots : angleMean (if not estimated)
#'
#' @return An object moveHMM.

moveHMM <- function(mod)
{
  if(is.null(mod$data) | is.null(mod$states) | is.null(mod$mle) | is.null(mod$stepDist) |
       is.null(mod$angleDist))
    stop("Can't construct moveHMM object : fields are missing")

  obj <- mod

  class(obj) <- append(class(obj),"moveHMM")
  return(obj)
}
