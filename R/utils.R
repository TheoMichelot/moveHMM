
#' Discrete colour palette for states
#'
#' @param nbStates Number of states
#'
#' @return Vector of colours, of length nbStates.
getPalette <- function(nbStates) {
    if(nbStates < 8) {
        # color-blind friendly palette
        pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        col <- pal[1:nbStates]
    } else {
        # to make sure that all colours are distinct (emulate ggplot default palette)
        hues <- seq(15, 375, length = nbStates + 1)
        col <- hcl(h = hues, l = 65, c = 100)[1:nbStates]
    }
    return(col)
}
