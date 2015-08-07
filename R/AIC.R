
#' AIC
#'
#' @param m A moveHMM object.
#'
#' @return The AIC of the fitted model.

AIC.moveHMM <- function(m)
{
  nbPar <- length(m$mle$stepPar)+length(m$mle$anglePar)+length(m$mle$beta)+length(m$mle$delta)-1
  maxLogLike <- -m$mod$minimum
  AIC <- 2*(nbPar-maxLogLike)
  return(AIC)
}
