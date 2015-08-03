#ifndef _DENSITIES_
#define _DENSITIES_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;

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

typedef arma::colvec (*FunPtr)(NumericVector,double,double);

#endif
