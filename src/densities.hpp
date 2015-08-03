#ifndef _DENSITIES_
#define _DENSITIES_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;

// Step length densities
arma::colvec dgamma_rcpp(NumericVector x, double mu, double sigma)
{
    arma::colvec res(x.size());
    double shape = pow(mu,2)/pow(sigma,2);
    double scale = pow(sigma,2)/mu;

    for(int i=0;i<x.size();i++) {
    	res(i) = R::dgamma(x(i),shape,scale,0);
	if(R_IsNA(res(i))) res(i) = 1;
    }

    return res;
}

arma::colvec dweibull_rcpp(NumericVector x, double mu, double sigma)
{
    arma::colvec res(x.size());
    double shape = pow(mu,2)/pow(sigma,2);
    double scale = pow(sigma,2)/mu;

    for(int i=0;i<x.size();i++) {
	res(i) = R::dweibull(x(i),shape,scale,0);
	if(R_IsNA(res(i))) res(i) = 1;
    }

    return res;
}

arma::colvec dexp_rcpp(NumericVector x, double rate, double foo=0)
{
     arma::colvec res(x.size());

    for(int i=0;i<x.size();i++) {
	res(i) = R::dexp(x(i),rate,0);
	if(R_IsNA(res(i))) res(i) = 1;
    }

    return res;
}

// Turning angle densities
arma::colvec dvm_rcpp(NumericVector x, double mu, double kappa)
{
    arma::colvec res(x.size());
    double b = R::bessel_i(kappa,0,2);

    for(int i=0;i<x.size();i++) {
	res(i) = 1/(2*M_PI*b)*pow((exp(cos(x(i)-mu)-1)),kappa);
	if(R_IsNA(res(i))) res(i) = 1;
    }

    return res;
}

arma::colvec dwrpcauchy_rcpp(NumericVector x, double mu, double rho)
{
    arma::colvec res(x.size());
    double num = 1-rho*rho;
    double den;

    for(int i=0;i<x.size();i++) {
	den = (2*M_PI)*(1+rho*rho-2*rho*cos(x(i)-mu));
	res(i) = num/den;
    }

    return res;
}

// used in nLogLike
typedef arma::colvec (*FunPtr)(NumericVector,double,double);

#endif
