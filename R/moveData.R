
#' Constructor of moveData objects
#'
#' @param data A list containing the IDs (if any), step lengths, turning angles
#' and covariates (if any), as well as the x/y values (either euclidean coordinates
#' or longitude/latitude)
#'
#' @return An object moveData.
#'
#' @examples
moveData <- function(data)
{
  if(is.null(data$step) | is.null(data$angle) | is.null(data$x) | is.null(data$y))
    stop("Can't construct moveData object : fields are missing")

  obj <- data

  class(obj) <- append(class(obj),"moveData")
  return(obj)
}
