
#' Summary \code{moveData}
#' @method summary moveData
#'
#' @param object A \code{moveData} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # Use example data set as output by prepData
#' summary(example$data)
#'
#' @export
summary.moveData <- function(object, ...)
{
    summary.data.frame(object, ...)
}
