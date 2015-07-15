
#' Constructor of moveData objects
#'
#' @param data A list containing the IDs (if any), step lengths, turning angles and covariates (if any).
#'
#' @return An object moveData.
#'
#' @examples
moveData <- function(data)
{
  obj <- data

  class(obj) <- append(class(obj),"moveData")
  return(obj)
}
