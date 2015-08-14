
#' Preprocessing of the tracking data
#'
#' @param trackData A dataframe of the tracking data, including at least the fields "x" and "y"
#' (either longitude/latitude values or cartesian coordinates), and optionnaly a field "ID"
#' (identifiers for the observed individuals). Additionnal fields are considered as covariates.
#' @param type 'GCD' if longitude/latitude provided (default), 'euclidean' otherwise.
#'
#' @return An object moveData, i.e. a dataframe of ID, step lengths,
#' turning angles and covariates values (if any).
#'
#' @examples
#' x <- c(1,2,3,4,5,6,7,8,9,10)
#' y <- c(1,1,1,2,2,2,1,1,1,2)
#' trackData <- data.frame(x=x,y=y)
#' d <- prepData(trackData,'euclidean')

prepData <- function(trackData, type=c('GCD','euclidean'))
{
  # check arguments
  type <- match.arg(type)
  x <- trackData$x
  y <- trackData$y
  if(is.null(x) | is.null(y))
    stop("The data does not have the right structure. x and y fields needed.")

  if(!is.null(trackData$ID)) ID <- as.character(trackData$ID) # homogenization of numeric and string IDs
  else ID <- rep("Animal1",length(x)) # default ID if none provided

  data <- data.frame(ID=character(),
                     step=numeric(),
                     angle=numeric())

  nbAnimals <- length(unique(ID))

  for(zoo in 1:nbAnimals) {
    nbObs <- length(which(ID==unique(ID)[zoo]))
    step <- rep(NA,nbObs)
    angle <- rep(NA,nbObs)
    i1 <- which(ID==unique(ID)[zoo])[1]
    i2 <- i1+nbObs-1
    for(i in (i1+1):(i2-1)) {
      if(!is.na(x[i-1]) & !is.na(x[i]) & !is.na(y[i-1]) & !is.na(y[i])) {
        step[i-i1+1] <- spDistsN1(pts = matrix(c(x[i-1],y[i-1]),ncol=2),
                                  pt = c(x[i],y[i]),
                                  longlat = (type=='GCD')) # TRUE if 'GCD', FALSE otherwise
      }

      if(!is.na(x[i-1]) & !is.na(x[i]) & !is.na(x[i+1]) & !is.na(y[i-1]) & !is.na(y[i]) & !is.na(y[i+1])) {
        angle[i-i1+1] <- turnAngle(c(x[i-1],y[i-1]),
                                   c(x[i],y[i]),
                                   c(x[i+1],y[i+1]))
      }
    }
    step[i2-i1+1] <- sqrt((x[i2]-x[i2-1])^2+(y[i2]-y[i2-1])^2)

    d <- data.frame(ID=rep(unique(ID)[zoo],nbObs),
                    step=step,
                    angle=angle)
    data <- rbind(data,d)
  }

  covsCol <- which(names(trackData)!="ID" & names(trackData)!="x" & names(trackData)!="y")
  if(length(covsCol)>0) {
    covs <- data.frame(trackData[,covsCol]) # to prevent error if nbCovs==1
    colnames(covs) <- names(trackData)[covsCol]

    # account for missing values of the covariates
    if(length(which(is.na(covs)))>0)
      warning(paste("There are",length(which(is.na(covs))),
                    "missing covariate values.",
                    "Each will be replaced by the closest available value."))
    for(i in 1:length(covsCol)) {
      if(length(which(is.na(covs[,i])))>0) { # if covariate i has missing values
        if(is.na(covs[1,i])) { # if the first value of the covariate is missing
          k <- 1
          while(is.na(covs[k,i])) k <- k+1
          for(j in k:2) covs[j-1,i] <- covs[j,i]
        }
        for(j in 2:nrow(trackData))
          if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
      }
    }
  }
  else covs <- NULL
  data <- cbind(data,x=trackData$x,y=trackData$y)
  if(!is.null(covs)) data <- cbind(data,covs)
  return(moveData(data))
}
