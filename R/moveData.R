
#' Constructor of moveData objects
#'
#' @param data A dataframe containing the IDs, step lengths, turning angles (if any),
#' and covariates (if any), as well as the x/y values (either cartesian coordinates
#' or longitude/latitude)
#'
#' @return An object moveData.

moveData <- function(data)
{
  if(is.null(data$step) | is.null(data$x) | is.null(data$y))
    stop("Can't construct moveData object : fields are missing")

  obj <- data

  class(obj) <- append(class(obj),"moveData")
  return(obj)
}
