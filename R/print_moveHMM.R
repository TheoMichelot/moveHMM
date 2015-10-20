
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
  p <- parDef(m$stepDist,m$angleDist,nbStates,TRUE,m$conditions$zeroInflation)

  if(length(m$mod)>1)
    cat("Value of the maximum log-likelihood:",-m$mod$minimum,"\n\n")

  cat("Step length parameters:\n")
  cat("----------------------\n")
  for(i in 1:nrow(m$mle$stepPar)) {
    cat(p$parNames[i],"\n")
    print(m$mle$stepPar[i,])
  }

  cat("\n")
  if(m$angleDist!="none") {
    cat("Turning angle parameters:\n")
    cat("------------------------\n")
    for(i in 1:nrow(m$mle$anglePar)) {
      cat(p$parNames[nrow(m$mle$stepPar)+i],"\n")
      print(m$mle$anglePar[i,])
    }
  }

  if(!is.null(m$mle$beta)) {
    cat("\n")
    cat("Regression coeffs for the transition probabilities:\n")
    cat("--------------------------------------------------\n")

    beta <- m$mle$beta
    f <- m$conditions$formula
    rownames(beta) <- c("intercept",attr(terms(f),"term.labels"))
    columns <- NULL
    for(i in 1:nbStates)
      for(j in 1:nbStates) {
        if(i<j)
          columns[(i-1)*nbStates+j-i] <- paste(i,"->",j)
        if(j<i)
          columns[(i-1)*(nbStates-1)+j] <- paste(i,"->",j)
      }
    colnames(beta) <- columns
    print(beta)
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
