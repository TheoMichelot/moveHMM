
#' Constructor of \code{moveData} objects
#'
#' This constructor is not intended to be used, except inside the function
#' \code{\link{prepData}}. Refer to the documentation for that function.
#'
#' @param data A data frame with columns: \code{ID} (track ID(s)), \code{step}
#' (step length), \code{angle} (turning angle, if any), \code{x} (Easting or
#' longitude), \code{y} (Norting or latitude), and covariates (if any).
#'
#' @return An object of class \code{moveData}.

moveData <- function(data)
{
    if(is.null(data$ID) | is.null(data$step) | is.null(data$x) | is.null(data$y))
        stop("Can't construct moveData object: fields are missing")

    obj <- data

    class(obj) <- append("moveData",class(obj))
    return(obj)
}

#' Is moveData
#'
#' Check that an object is of class \code{\link{moveData}}. Used in \code{\link{fitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{moveData}}, \code{FALSE} otherwise.

is.moveData <- function(x)
    inherits(x,"moveData")
