
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

#' Density function of von Mises distribution
#'
#' @param x Angle
#' @param mu Mean parameter
#' @param kappa Concentration parameter
#' @param log Should log-density be returned?
#'
#' @return Von Mises density
dvm <- function(x, mu, kappa, log = FALSE) {
    # The "- kappa" term below cancels out the expon.scaled
    b <- besselI(kappa, 0, expon.scaled = TRUE)
    val <- - log(2 * pi * b) + kappa * cos(x - mu) - kappa
    if(!log) {
        val <- exp(val)
    }
    return(val)
}

#' Sample from von Mises distribution
#'
#' @param n Number of samples
#' @param mu Mean parameter
#' @param kappa Concentration parameter
#'
#' @return Vector of n samples from vm(mu, kappa)
#'
#' @details Uses basic rejection sampling, based on dvm(), which might
#' be inefficient for large kappa. Could be improved following Best & Fisher
#' (1979), Efficient simulation of the von Mises distribution, JRSSC, 28(2),
#' 152-157.
rvm <- function(n, mu, kappa) {
    x <- rep(NA, n)
    n_accept <- 0
    pdf_max <- dvm(x = mu, mu = mu, kappa = kappa)
    while(n_accept < n) {
        x_star <- runif(1, min = -pi, max = pi)
        pdf_star <- dvm(x = x_star, mu = mu, kappa = kappa)
        accept_prob <- pdf_star/pdf_max
        if(runif(1) < accept_prob) {
            n_accept <- n_accept + 1
            x[n_accept] <- x_star
        }
    }
    return(x)
}

#' Density function of wrapped Cauchy distribution
#'
#' @param x Angle
#' @param mu Mean parameter
#' @param rho Concentration parameter
#' @param log Should log-density be returned?
#'
#' @return Wrapped Cauchy density
dwrpcauchy <- function(x, mu, rho, log = FALSE) {
    val <- (1 - rho^2) /
        (2 * pi * (1 + rho^2 - 2 * rho * cos(x - mu)))
    if(log) {
        val <- log(val)
    }
    return(val)
}

#' Sample from wrapped Cauchy distribution
#'
#' @param n Number of samples
#' @param mu Mean parameter
#' @param rho Concentration parameter
#'
#' @return Vector of n samples from wrpcauchy(mu, rho)
#'
#' @details Uses basic rejection sampling, based on dwrpcauchy(), which might
#' be inefficient for large rho.
rwrpcauchy <- function(n, mu, rho) {
    x <- rep(NA, n)
    n_accept <- 0
    pdf_max <- dwrpcauchy(x = mu, mu = mu, rho = rho)
    while(n_accept < n) {
        x_star <- runif(1, min = -pi, max = pi)
        pdf_star <- dwrpcauchy(x = x_star, mu = mu, rho = rho)
        accept_prob <- pdf_star/pdf_max
        if(runif(1) < accept_prob) {
            n_accept <- n_accept + 1
            x[n_accept] <- x_star
        }
    }
    return(x)
}
