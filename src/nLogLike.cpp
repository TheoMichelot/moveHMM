#include "densities.h"

// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat beta, arma::mat covs, DataFrame data, std::string stepDist, 
                        std::string angleDist, arma::mat stepPar, arma::mat anglePar, 
                        arma::rowvec delta, IntegerVector aInd, bool zeroInflation=false,
                        bool stationary=false)
{
    // 1. Computation of transition probability matrix trMat
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
    
    // 2. Computation of matrix of joint probabilities allProbs
    
    // map the functions names with the actual functions
    // (the type FunPtr and the density functions are defined in densities.h)
    map<std::string,FunPtr> funMap;
    funMap["gamma"] = dgamma_rcpp;
    funMap["weibull"] = dweibull_rcpp;
    funMap["lnorm"] = dlnorm_rcpp;
    funMap["exp"] = dexp_rcpp;
    funMap["vm"] = dvm_rcpp;
    funMap["wrpcauchy"] = dwrpcauchy_rcpp;

    // compute stationary distribution delta
    if(stationary) { 
        arma::mat diag(nbStates,nbStates);
        diag.eye(); // diagonal of ones
        arma::mat gamma = trMat.slice(0); // all slices are identical if stationary
        arma::colvec v(nbStates);
        v.ones(); // vector of ones
        delta = arma::solve(diag-gamma+1,v).t();
    }
    
    arma::mat allProbs(nbObs,nbStates);
    allProbs.ones();
    
    arma::colvec stepProb(nbObs);
    NumericVector stepArgs(2);
    NumericVector step = data["step"];
   
    arma::colvec angleProb(nbObs);
    NumericVector angleArgs(2);
    NumericVector angle(nbObs);
    
    if(angleDist!="none") {
        angle = data["angle"];
    }

    arma::rowvec zeromass(nbStates);

    if(zeroInflation) {
	    zeromass = stepPar.row(stepPar.n_rows-1);
	    arma::mat stepPar2 = stepPar.submat(0,0,stepPar.n_rows-2,stepPar.n_cols-1);
	    stepPar = stepPar2;
    }

    for(int state=0;state<nbStates;state++) 
    {
        for(int i=0;i<stepPar.n_rows;i++)
            stepArgs(i) = stepPar(i,state);

        if(zeroInflation) {
            // remove the NAs from step (impossible to subset a vector with NAs)
            for(int i=0;i<nbObs;i++) {
                if(R_IsNA(step(i))) {
                    step(i) = -1;
                    stepProb(i) = 1;
                }
            }

            // compute probability of non-zero observations
            stepProb.elem(arma::find(as<arma::vec>(step)>0)) = (1-zeromass(state))*funMap[stepDist](step[step>0],stepArgs(0),stepArgs(1));
            
            // compute probability of zero observations
            int nbZeros = as<NumericVector>(step[step==0]).size();
            arma::vec zm(nbZeros);
            for(int i=0;i<nbZeros;i++)
                zm(i) = zeromass(state);
            stepProb.elem(arma::find(as<arma::vec>(step)==0)) = zm;

            // put the NAs back
            for(int i=0;i<nbObs;i++) {
                if(step(i)<0)
                    step(i) = NA_REAL;
            }
        }
        else
	        stepProb = funMap[stepDist](step,stepArgs(0),stepArgs(1));

	    if(angleDist!="none") {
	        for(int i=0;i<anglePar.n_rows;i++)
		        angleArgs(i) = anglePar(i,state);

	        angleProb = funMap[angleDist](angle,angleArgs(0),angleArgs(1));
	        allProbs.col(state) = stepProb%angleProb;
	    }
	    else allProbs.col(state) = stepProb;
    }

    for(int i=0;i<30;i++)
        cout << allProbs(i,0) << ", " << allProbs(i,1) << endl;

    // 3. Forward algorithm
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

