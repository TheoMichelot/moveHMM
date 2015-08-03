#include "densities.hpp"

// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat beta, arma::mat covs, DataFrame data, std::string stepDist, 
                        std::string angleDist, arma::mat stepPar, arma::mat anglePar, 
                        arma::rowvec delta, IntegerVector aInd, bool zeroInflation=false)
{
    // Computation of transition probability matrix trMat
    int nbObs = data.nrows();
    arma::cube trMat(nbStates,nbStates,nbObs);
    trMat.zeros();
    arma::mat rowSums(nbStates,nbObs);
    rowSums.zeros();

    arma::mat g(nbObs,nbStates*(nbStates-1));
    g = covs*beta;

    for(int k=0;k<nbObs;k++) {
        int cpt=0;
        for(int i=0;i<nbStates;i++) {
            for(int j=0;j<nbStates;j++) {
                if(i==j) {
                    trMat(i,j,k)=1;
                    cpt++;
                }
                else trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));
                rowSums(i,k)=rowSums(i,k)+trMat(i,j,k);
            }
        }
    }

    // normalization
    for(int k=0;k<nbObs;k++)
        for(int i=0;i<nbStates;i++)
            for(int j=0;j<nbStates;j++)
                trMat(i,j,k) = trMat(i,j,k)/rowSums(i,k);
    
    // Computation of matrix of joint probabilities allProbs
    map<std::string,FunPtr> funMap;
    funMap["gamma"]=dgamma_rcpp;
    funMap["weibull"]=dweibull_rcpp;
    funMap["exp"]=dexp_rcpp;
    funMap["vm"]=dvm_rcpp;
    funMap["wrpcauchy"]=dwrpcauchy_rcpp;
	
    arma::mat allProbs(nbObs,nbStates);
    allProbs.ones();
    
    arma::colvec stepProb(nbObs);
    NumericVector stepArgs(2);
    NumericVector step = data["step"];
   
    arma::colvec angleProb(nbObs);
    NumericVector angleArgs(2);
    NumericVector angle(nbObs);

    if(angleDist!="NULL")
	angle = data["angle"];

    arma::rowvec zeromass(nbStates);

    if(zeroInflation) {
	zeromass = stepPar.row(stepPar.n_rows-1);
	arma::mat stepPar = stepPar(arma::span(0,stepPar.n_rows-2),arma::span());
    }

    for(int state=0;state<nbStates;state++) 
    {
        for(int i=0;i<stepPar.n_rows;i++)
            stepArgs(i) = stepPar(i,state);
	
	stepProb = funMap[stepDist](step,stepArgs(0),stepArgs(1));
	if(zeroInflation) {
	    for(int i=0;i<nbObs;i++) {
		if(step(i)==0) stepProb(i)=zeromass(state);
	    }
	}

	if(angleDist!="NULL") {
	    for(int i=0;i<anglePar.n_rows;i++)
		angleArgs(i) = anglePar(i,state);

	    angleProb = funMap[angleDist](angle,angleArgs(0),angleArgs(1));
	    allProbs.col(state) = stepProb%angleProb;
	}
	else allProbs.col(state) = stepProb;
    
    }

    // Forward algorithm
    arma::mat gamma(nbStates,nbStates);
    double lscale = 0;
    int k=1;
    arma::rowvec alpha = delta%allProbs.row(0);

    for(int i=1;i<allProbs.n_rows;i++) {
	if(k<aInd.size() && i==aInd(k)-1) {
	    k++;
	    alpha = delta%allProbs.row(i);
	}
        gamma = trMat.slice(i);
        alpha = alpha*gamma%allProbs.row(i);

        lscale = lscale + log(sum(alpha));
        alpha = alpha/sum(alpha);
    }

    return -lscale;
}

