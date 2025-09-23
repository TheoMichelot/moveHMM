#ifndef _DENSITIES_
#define _DENSITIES_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;

//' Gamma density function
//'
//' Probability density function of the gamma distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param sigma Standard deviation
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dgamma_rcpp(NumericVector x, double mu, double sigma)
{
    arma::colvec res(x.size());

    // convert mean and sd to shape and scale
    double shape = pow(mu,2)/pow(sigma,2);
    double scale = pow(sigma,2)/mu;

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else
            res(i) = R::dgamma(x(i),shape,scale,0);
    }

    return res;
}

//' Weibull density function
//'
//' Probability density function of the Weibull distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param shape Shape
//' @param scale Scale
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dweibull_rcpp(NumericVector x, double shape, double scale)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else
            res(i) = R::dweibull(x(i),shape,scale,0);
    }

    return res;
}

//' Log-normal density function
//'
//' Probability density function of the log-normal distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param meanlog Mean of the distribution on the log-scale
//' @param sdlog Standard deviation of the distribution on the log-scale
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dlnorm_rcpp(NumericVector x, double meanlog, double sdlog)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else
            res(i) = R::dlnorm(x(i),meanlog,sdlog,0);
    }

    return res;
}

//' Exponential density function
//'
//' Probability density function of the exponential distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param rate Rate
//' @param foo Unused (for compatibility with template)
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dexp_rcpp(NumericVector x, double rate, double foo=0)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else
            res(i) = R::dexp(x(i),1/rate,0); // R::dexp expects scale=1/rate
    }

    return res;
}

//' Von Mises density function
//'
//' Probability density function of the Von Mises distribution, defined as a function
//' of the modified Bessel function of order 0 (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param kappa Concentration
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dvm_rcpp(NumericVector x, double mu, double kappa)
{
    arma::colvec res(x.size());
    double b = R::bessel_i(kappa, 0, 2);

    // In the formula below, there is an additional exp(-kappa) factor
    // (compared to the standard von Mises pdf), because b is calculated
    // with the option expon.scaled = TRUE, i.e. multiplied by exp(-kappa).
    // See the documentation for besselI in R for details.

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else
            res(i) = 1/(2 * M_PI * b) * exp(kappa * cos(x(i) - mu) - kappa);
    }

    return res;
}

//' Wrapped Cauchy density function
//'
//' Probability density function of the wrapped Cauchy distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param rho Concentration
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dwrpcauchy_rcpp(NumericVector x, double mu, double rho)
{
    arma::colvec res(x.size());
    double num = 1-rho*rho;
    double den;

    for(int i=0;i<x.size();i++) {
        if(!std::isfinite(x(i)))
            res(i) = 1; // if missing observation
        else {
            den = (2*M_PI)*(1+rho*rho-2*rho*cos(x(i)-mu));
            res(i) = num/den;
        }
    }

    return res;
}

// used in nLogLike_rcpp to map the functions' names to the functions
typedef arma::colvec (*FunPtr)(NumericVector,double,double);

#endif
