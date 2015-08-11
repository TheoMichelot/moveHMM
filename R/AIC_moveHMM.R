
#' AIC
#'
#' @param object A moveHMM object.
#' @param ... Currently unused. For compatibility with generic method.
#' @param k Penalty per parameter (default k=2 for classical AIC).
#'
#' @return The AIC of the fitted model.

AIC.moveHMM <- function(object,...,k=2)
{
  m <- object
  nbPar <- length(m$mle$stepPar)+length(m$mle$anglePar)+length(m$mle$beta)+length(m$mle$delta)-1
  maxLogLike <- -m$mod$minimum
  AIC <- -2*maxLogLike+k*nbPar
  return(AIC)
}
