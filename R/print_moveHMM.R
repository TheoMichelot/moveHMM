
#' Print moveHMM
#' @method print moveHMM
#'
#' @param x A moveHMM object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#'
#' print(m)

print.moveHMM <- function(x,...)
{
  m <- x
  nbStates <- ncol(m$mle$stepPar)
  p <- parDef(m$stepDist,m$angleDist,nbStates,TRUE,m$conditions$zeroInflation)

  cat("Value of the maximum log-likelihood :",-m$mod$minimum,"\n\n")

  cat("Step length parameters :\n")
  cat("----------------------\n")
  for(i in 1:nrow(m$mle$stepPar)) {
    cat(p$parNames[i],"\n")
    print(m$mle$stepPar[i,])
  }

  cat("\n")
  if(m$angleDist!="none") {
    cat("Turning angle parameters :\n")
    cat("------------------------\n")
    for(i in 1:nrow(m$mle$anglePar)) {
      cat(p$parNames[nrow(m$mle$stepPar)+i],"\n")
      print(m$mle$anglePar[i,])
    }
  }

  cat("\n")
  cat("Transition probabilities parameters :\n")
  cat("-----------------------------------\n")
  print(m$mle$beta)

  cat("\n")
  cat("Initial distribution :\n")
  cat("--------------------\n")
  print(m$mle$delta)
}
