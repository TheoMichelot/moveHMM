
#' Turning angle
#'
#' Used in \code{\link{prepData}}.
#'
#' @param x First point
#' @param y Second point
#' @param z Third point
#' @param LLangle Logical. If TRUE, the turning angle is calculated with
#' \code{geosphere::bearing}, else calculated with \code{atan2}.
#'
#' @return The angle between vectors (x,y) and (y,z)
#'
#' @examples
#' \dontrun{
#' x <- c(0,0)
#' y <- c(4,6)
#' z <- c(10,7)
#' turnAngle(x,y,z,LLangle=FALSE)
#' }
#'
#' @importFrom geosphere bearing

turnAngle <- function(x,y,z,LLangle)
{
    # NA angle if zero step length
    if(all(x==y) | all(y==z))
        return(NA)

    if(LLangle) {
        angle <- (bearing(x,y)-bearing(y,z))/180*pi
    } else {
        v <- c(y[1]-x[1],y[2]-x[2])
        w <- c(z[1]-y[1],z[2]-y[2])
        angle <- atan2(w[2],w[1])-atan2(v[2],v[1])
    }

    while(angle<=(-pi))
        angle <- angle + 2*pi
    while(angle>pi)
        angle <- angle - 2*pi

    return(angle)
}
