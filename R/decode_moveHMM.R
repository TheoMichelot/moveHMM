
#' Generic decode method
#' @param m Fitted model
decode <- function(m) UseMethod("decode") # define generic method decode

#' Decode the model
#'
#' @param m moveHMM object
#'
#' @return List of : sp (state probabilities), ...
#' @examples
#' m <- ex$mod # moveHMM object (returned by fitHMM)
#'
#' dec <- decode(m)

decode.moveHMM <- function(m)
{
  # states probabilities
  sp <- stateProbs(m)

  # pseudo-residuals
  res <- pseudoRes(m)

  d <- list(sp=sp,res=res)
  return(d)
}
