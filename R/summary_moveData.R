
#' Summary \code{moveData}
#' @method summary moveData
#'
#' @param object A \code{moveData} object.
#' @param details \code{TRUE} if quantiles of the covariate values should be printed,
#' \code{FALSE} otherwise (default).
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' data <- ex$data # moveData object (returned by prepData or simData)
#'
#' summary(data)
#'
#' @export

summary.moveData <- function(object,details=FALSE,...)
{
  data <- object

  # print animals' IDs and numbers of observations
  nbAnimals <- length(unique(data$ID))
  if(nbAnimals==1)
    cat("Movement data for 1 animal:\n",sep="")
  else
    cat("Movement data for ",nbAnimals," animals:\n",sep="")
  for(zoo in 1:nbAnimals)
    cat(as.character(unique(data$ID)[zoo])," -- ",
        length(which(data$ID==unique(data$ID)[zoo]))," observations\n",sep="")

  # identify covariates
  covsCol <- which(names(data)!="ID" & names(data)!="x" & names(data)!="y" &
                     names(data)!="step" & names(data)!="angle" & names(data)!="states")

  # print covariates' names
  if(length(covsCol)>0) {
    cat("\nCovariate(s): ",sep="")
    if(!details) {
      if(length(covsCol)>1) {
        for(i in 1:(length(covsCol)-1))
          cat(names(data)[covsCol[i]],", ",sep="")
      }
      cat(names(data)[covsCol[length(covsCol)]])
    } else {
      # print names and quantiles of the covariates
      for(i in 1:length(covsCol)) {
        cat("\n",names(data)[covsCol[i]],"\n")
        q <- quantile(data[,covsCol[i]])
        # insert the mean in the quantiles' vector
        q[6] <- q[5]
        q[5] <- q[4]
        q[4] <- mean(data[,covsCol[i]])
        names(q) <- c("Min.","25%","Median","Mean","75%","Max.")
        print(q)
      }
    }
  } else {
    cat("No covariates.\n")
  }
}
