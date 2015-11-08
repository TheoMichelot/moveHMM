
#' Print \code{moveHMM}
#' @method print moveHMM
#'
#' @param x A \code{moveHMM} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' # m is a moveHMM object (as returned by fitHMM), automatically loaded with the package
#' m <- example$m
#'
#' print(m)
#'
#' @export

print.moveHMM <- function(x,...)
{
  m <- x
  nbStates <- ncol(m$mle$stepPar)
  p <- parDef(m$conditions$stepDist,m$conditions$angleDist,nbStates,TRUE,
              m$conditions$zeroInflation)

  if(length(m$mod)>1)
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")

  cat("Step length parameters:\n")
  cat("----------------------\n")
  print(m$mle$stepPar)

  cat("\n")
  if(m$conditions$angleDist!="none") {
    cat("Turning angle parameters:\n")
    cat("------------------------\n")
    print(m$mle$anglePar)
  }

  if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("--------------------------------------------------\n")
    print(m$mle$beta)
  }

  if(!is.null(m$mle$gamma)) {
    cat("\n")
    cat("Transition probability matrix:\n")
    cat("-----------------------------\n")
    print(m$mle$gamma)
  }

  cat("\n")
  cat("Initial distribution:\n")
  cat("--------------------\n")
  print(m$mle$delta)
}
