
#' Print moveHMM
#'
#' @param m A moveHMM object.
#'
#' @examples
#' m <- example$mod # moveHMM object (returned by fitHMM)
#'
#' print(m)

print.moveHMM <- function(m)
{
  nbStates <- ncol(m$mle$stepPar)
  p <- parDef(m$stepDist,m$angleDist,nbStates,TRUE,m$zeroInflation)

  cat("Step length parameters :\n")
  cat("----------------------\n")
  for(i in 1:nrow(m$mle$stepPar))
    cat(p$parNames[i],"\t",m$mle$stepPar[i,],"\n")

  cat("\n")
  if(angleDist!="none") {
    cat("Turning angle parameters :\n")
    cat("------------------------\n")
    for(i in 1:nrow(m$mle$anglePar))
      cat(p$parNames[nrow(m$mle$stepPar)+i],"\t",m$mle$anglePar[i,],"\n")
  }

  cat("\n")
  cat("Beta :\n")
  cat("----\n")
  print(m$mle$beta)

  cat("\n")
  cat("Initial distribution :\n")
  cat("--------------------\n")
  cat(m$mle$delta,"\n")
}
