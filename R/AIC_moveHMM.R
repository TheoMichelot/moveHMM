
#' AIC
#'
#' Akaike information criterion of a moveHMM model.
#'
#' @method AIC moveHMM
#'
#' @param object A moveHMM object.
#' @param ... Optional additional moveHMM objects, to compare AICs of the different models.
#' @param k Penalty per parameter. Default : 2 ; for classical AIC.
#'
#' @return The AIC of the model(s) provided.

AIC.moveHMM <- function(object,...,k=2)
{
  models <- list(...)

  if(length(models)>0) { # if several models are provided
    modNames <- all.vars(match.call()) # store the names of the models given as arguments

    models[[length(models)+1]] <- object
    AIC <- rep(NA,length(models))

    for(i in 1:length(models)) {
      m <- models[[i]]
      nbPar <- length(m$mle$stepPar)+length(m$mle$anglePar)+length(m$mle$beta)+length(m$mle$delta)-1
      maxLogLike <- -m$mod$minimum
      AIC[i] <- -2*maxLogLike+k*nbPar
    }

    return(data.frame(Model=modNames,AIC=AIC))
  }
  else { # if only one model is provided
    m <- object
    nbPar <- length(m$mle$stepPar)+length(m$mle$anglePar)+length(m$mle$beta)+length(m$mle$delta)-1
    maxLogLike <- -m$mod$minimum
    AIC <- -2*maxLogLike+k*nbPar

    return(AIC)
  }
}
