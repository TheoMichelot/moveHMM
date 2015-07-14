
#' Constructor of moveData objects
#'
#' @param ID Vector of the IDs of observed individuals
#' @param stepLength Vector of step lengths
#' @param turnAngle Vector of turning angles
#' @param covs Matrix of covariates if any
#'
#' @return An object moveData, i.e. a list containing the IDs, the lengths, the angles and the covariates.
#'
#' @examples
moveData <- function(ID,length,angle,covs=NULL)
{
  obj <- list(
    ID = ID,
    stepLength = length,
    turnAngle = angle,
    covs = covs
  )

  class(obj) <- append(class(obj),"moveData")
  return(obj)
}
