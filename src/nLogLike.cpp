#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double nLogLike_rcpp(NumericVector trVec, arma::rowvec delta, arma::mat allProbs)
{
    IntegerVector dims = trVec.attr("dim");
    arma::cube trMat(trVec.begin(),dims[0],dims[1],dims[2],false);
    int nbStates = trMat.slice(0).n_rows;
    arma::mat gamma(nbStates,nbStates);
    double lscale = 0;
    arma::rowvec alpha = delta%allProbs.row(0);

    for(int i=1;i<allProbs.n_rows;i++) {
	gamma = trMat.slice(i);
	alpha = alpha*gamma%allProbs.row(i);

	lscale = lscale + log(sum(alpha));
	alpha = alpha/sum(alpha);
    }

    return lscale;
}

