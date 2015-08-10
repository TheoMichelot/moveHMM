
decode <- function(m) UseMethod("decode") # define generic method decode

#' Decode the model
#'
#' @param m moveHMM object
#'
#' @return List of : sp (state probabilities), ...
#' @examples
#' m <- example$mod # moveHMM object (returned by fitHMM)
#'
#' dec <- decode(m)

decode.moveHMM <- function(m)
{
  # states probabilities
  sp <- stateProbs(m)

  #residuals...
  d <- list(sp=sp)
  return(d)
}
