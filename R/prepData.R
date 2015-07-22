
#' Preprocessing of the tracking data
#'
#' @param trackData A dataframe of the tracking data, including at least the fields "x" and "y"
#' (either longitude/latitude values or euclidean coordinates), and optionnaly a field "ID" (identifiers
#' for the observed individuals). Additionnal fields are considered as covariates.
#' @param type 'GCD' if longitude/latitude provided (default), 'euclidean' otherwise.
#'
#' @return An object moveData, i.e. a list in which each element is a list of ID, step lengths,
#' turning angles and covariates values corresponding to an individual.
#'
#' @examples
#' x <- c(1,2,3,4,5,6,7,8,9,10)
#' y <- c(1,1,1,2,2,2,1,1,1,2)
#' trackData <- list(x,y)
#' d <- prepData(trackData,'euclidean')
#' d[[1]]$step # the vector of step lengths of the first individual

prepData <- function(trackData, type=c('GCD','euclidean'))
{
  type <- match.arg(type)
  x <- trackData$x
  y <- trackData$y
  if(!is.null(trackData$ID)) ID <- as.character(trackData$ID) # homogenization of numeric and string IDs
  else ID <- NULL
  covsCol <- which(names(trackData)!="ID" & names(trackData)!="x" & names(trackData)!="y")
  if(length(covsCol)>0) covs <- trackData[,covsCol]
  else covs <- NULL

  if(is.null(x) | is.null(y)) stop("The data does not have the right structure.")

  data <- list() # to be returned

  if(is.null(ID)) nbAnimals <-1
  else nbAnimals <- length(unique(ID))

  for(k in 1:nbAnimals) {
    if(is.null(ID)) nbObs <- length(x)
    else nbObs <- length(which(ID==unique(ID)[k]))
    step <- rep(NA,nbObs)
    angle <- rep(NA,nbObs)
    if(!is.null(ID)) i1 <- which(ID==unique(ID)[k])[1]
    else i1 <- 1
    i2 <- i1+nbObs-1
    for(i in (i1+1):(i2-1)) {
      step[i-i1+1] <- spDistsN1(pts = matrix(c(x[i-1],y[i-1]),ncol=2),
                              pt = c(x[i],y[i]),
                              longlat = (type=='GCD')) # TRUE if 'GCD', FALSE otherwise

      angle[i-i1+1] <- turnAngle(c(x[i-1],y[i-1]),
                                 c(x[i],y[i]),
                                 c(x[i+1],y[i+1]))
    }
    step[i2-i1+1] <- sqrt((x[i2]-x[i2-1])^2+(y[i2]-y[i2-1])^2)

    c <- NULL
    if(!is.null(covs)) {
      if(length(covsCol)==1) c <- covs[i1:i2]
      else c <- covs[i1:i2,]
    }

    data[[k]] <- list(ID=unique(ID)[k],step=step,angle=angle,covs=c,x=x[i1:i2],y=y[i1:i2])
  }

  return(moveData(data))
}
