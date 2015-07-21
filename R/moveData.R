
#' Constructor of moveData objects
#'
#' @param data A list containing the IDs (if any), step lengths, turning angles
#' and covariates (if any), as well as the x/y values (either euclidean coordinates
#' or longitude/latitude)
#'
#' @return An object moveData.

moveData <- function(data)
{
  for(i in 1:length(data))
    if(is.null(data[[i]]$step) | is.null(data[[i]]$angle) | is.null(data[[i]]$x) | is.null(data[[i]]$y))
      stop("Can't construct moveData object : fields are missing")

  obj <- data

  class(obj) <- append(class(obj),"moveData")
  return(obj)
}
