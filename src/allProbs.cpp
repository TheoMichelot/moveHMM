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

// [[Rcpp::export]]
arma::mat allProbs_rcpp(DataFrame data, int nbStates, std::string stepDist, std::string angleDist, arma::mat stepPar, 
						arma::mat anglePar, bool zeroInflation=false)
{
    map<std::string,FunPtr> funMap;
    funMap["gamma"]=dgamma_rcpp;
    funMap["vm"]=dvm_rcpp;
	
    int nbObs = data.nrows();
    arma::mat allProbs(nbObs,nbStates);
    allProbs.ones();

    arma::colvec stepProb(nbObs);
    arma::colvec angleProb(nbObs);
    NumericVector stepArgs(2);
    NumericVector angleArgs(2);

    for(int state=0;state<nbStates;state++) 
    {
	for(int i=0;i<stepPar.n_rows;i++)
	    stepArgs(i) = stepPar(i,state);
	for(int i=0;i<anglePar.n_rows;i++)
	    angleArgs(i) = anglePar(i,state);

	stepProb = funMap[stepDist](data["step"],stepArgs(0),stepArgs(1));
	angleProb = funMap[angleDist](data["angle"],angleArgs(0),angleArgs(1));

	allProbs.col(state) = stepProb%angleProb;
    }

    return allProbs;
}
