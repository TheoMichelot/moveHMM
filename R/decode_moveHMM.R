
#' Generic decode method
#' @param m Fitted model
decode <- function(m) UseMethod("decode") # define generic method decode

#' Decode the model
#' @method decode moveHMM
#'
#' @param m moveHMM object
#'
#' @return List of : sp (state probabilities), res (pseudo-residuals),...
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#'
#' dec <- decode(m)
#' dec$sp # state probabilities
#' dec$res # pseudo-residuals

decode.moveHMM <- function(m)
{
  # states probabilities
  sp <- stateProbs(m)

  # pseudo-residuals
  res <- pseudoRes(m)

  d <- list(sp=sp,res=res)
  return(d)
}
