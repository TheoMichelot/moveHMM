#ifndef _DENSITIES_
#define _DENSITIES_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;

// Probability density function of the Gamma distribution
// [[Rcpp::export]]
arma::colvec dgamma_rcpp(NumericVector x, double mu, double sigma)
{
    arma::colvec res(x.size());

    // convert mean and sd to shape and scale
    double shape = pow(mu,2)/pow(sigma,2);
    double scale = pow(sigma,2)/mu;

    for(int i=0;i<x.size();i++) {
    	res(i) = R::dgamma(x(i),shape,scale,0);
	    if(R_IsNA(res(i))) 
            res(i) = 1; // if missing observation
    }

    return res;
}

// Probability density function of the Weibull distribution
// [[Rcpp::export]]
arma::colvec dweibull_rcpp(NumericVector x, double shape, double scale)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
	    res(i) = R::dweibull(x(i),shape,scale,0);
	    if(R_IsNA(res(i))) 
            res(i) = 1; // if missing observation
    }

    return res;
}

// Probability density function of the log-normal distribution
// [[Rcpp::export]]
arma::colvec dlnorm_rcpp(NumericVector x, double meanlog, double sdlog)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
	    res(i) = R::dlnorm(x(i),meanlog,sdlog,0);
	    if(R_IsNA(res(i))) 
            res(i) = 1; // if missing observation
    }

    return res;
}

// Probability density function of the exponential distribution
// [[Rcpp::export]]
arma::colvec dexp_rcpp(NumericVector x, double rate, double foo=0)
{
    arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
	    res(i) = R::dexp(x(i),1/rate,0); // R::dexp expects scale=1/rate
	    if(R_IsNA(res(i))) 
            res(i) = 1; // if missing observation
    }

    return res;
}

// Probability density function of the Von Mises distribution
// (defined as a function of the modified Bessel function of order 0
// [[Rcpp::export]]
arma::colvec dvm_rcpp(NumericVector x, double mu, double kappa)
{
    arma::colvec res(x.size());
    double b = R::bessel_i(kappa,0,2);

    for(int i=0;i<x.size();i++) {
	    res(i) = 1/(2*M_PI*b)*pow((exp(cos(x(i)-mu)-1)),kappa);
	    if(R_IsNA(res(i))) 
            res(i) = 1; // is missing observation
    }

    return res;
}

// Probability density function of the wrapped Cauchy distribution
// [[Rcpp::export]]
arma::colvec dwrpcauchy_rcpp(NumericVector x, double mu, double rho)
{
    arma::colvec res(x.size());
    double num = 1-rho*rho;
    double den;

    for(int i=0;i<x.size();i++) {
	    den = (2*M_PI)*(1+rho*rho-2*rho*cos(x(i)-mu));
	    res(i) = num/den;
	    if(R_IsNA(res(i))) 
            res(i) = 1; // if missing observation
    }

    return res;
}

// used in nLogLike_rcpp to map the functions' names to the functions
typedef arma::colvec (*FunPtr)(NumericVector,double,double);

#endif
