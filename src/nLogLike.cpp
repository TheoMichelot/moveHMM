#include "densities.h"

//' Negative log-likelihood
//'
//' Computation of the negative log-likelihood (forward algorithm - written in C++)
//'
//' @param nbStates Number of states
//' @param beta Matrix of regression coefficients for the transition probabilities
//' @param covs Covariates
//' @param data A \code{\link{moveData}} object of the observations
//' @param stepDist The name of the step length distribution
//' @param angleDist The name of the turning angle distribution
//' @param stepPar State-dependent parameters of the step length distribution
//' @param anglePar State-dependent parameters of the turning angle distribution
//' @param delta Stationary distribution
//' @param aInd Vector of indices of the rows at which the data switches to another animal
//' @param zeroInflation \code{true} if zero-inflation is included in the step length distribution,
//' \code{false} otherwise.
//' @param stationary \code{false} if there are covariates. If \code{true}, the initial distribution
//' is considered equal to the stationary distribution.
//' @param knownStates Vector of values of the state process which are known prior to fitting the
//' model (if any). Default: NULL (states are not known). This should be a vector with length the number
//' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
//' the state is not known.
//'
//' @return Negative log-likelihood
// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat beta, arma::mat covs, DataFrame data, std::string stepDist,
                     std::string angleDist, arma::mat stepPar, arma::mat anglePar, arma::rowvec delta,
                     IntegerVector aInd, bool zeroInflation, bool stationary, IntegerVector knownStates)
{
    int nbObs = data.nrows();

    //=======================================================//
    // 1. Computation of transition probability matrix trMat //
    //=======================================================//

    arma::cube trMat(nbStates,nbStates,nbObs);
    trMat.zeros();
    arma::mat rowSums(nbStates,nbObs);
    rowSums.zeros();

    arma::mat g(nbObs,nbStates*(nbStates-1));

    if(nbStates>1) {
        g = covs*beta;

        for(int k=0;k<nbObs;k++) {
            int cpt=0; // counter for diagonal elements
            for(int i=0;i<nbStates;i++) {
                for(int j=0;j<nbStates;j++) {
                    if(i==j) {
                        // if diagonal element, set to one and increment counter
                        trMat(i,j,k)=1;
                        cpt++;
                    }
                    else
                        trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));

                    // keep track of row sums, to normalize in the end
                    rowSums(i,k)=rowSums(i,k)+trMat(i,j,k);
                }
            }
        }

        // normalization
        for(int k=0;k<nbObs;k++)
            for(int i=0;i<nbStates;i++)
                for(int j=0;j<nbStates;j++)
                    trMat(i,j,k) = trMat(i,j,k)/rowSums(i,k);
    }

    //==========================================================//
    // 2. Computation of matrix of joint probabilities allProbs //
    //==========================================================//

    // map the functions names with the actual functions
    // (the type FunPtr and the density functions are defined in densities.h)
    map<std::string,FunPtr> funMap;
    funMap["gamma"] = dgamma_rcpp;
    funMap["weibull"] = dweibull_rcpp;
    funMap["lnorm"] = dlnorm_rcpp;
    funMap["exp"] = dexp_rcpp;
    funMap["vm"] = dvm_rcpp;
    funMap["wrpcauchy"] = dwrpcauchy_rcpp;

    if(nbStates==1)
        delta = 1; // no distribution if only one state
    else if(stationary) {
        // compute stationary distribution delta

        arma::mat diag(nbStates,nbStates);
        diag.eye(); // diagonal of ones
        arma::mat Gamma = trMat.slice(0).t(); // all slices are identical if stationary
        arma::colvec v(nbStates);
        v.ones(); // vector of ones
        try {
            delta = arma::solve(diag-Gamma+1,v).t();
        }
        catch(...) {
            throw std::runtime_error("A problem occurred in the calculation of "
                                         "the stationary distribution. You may want to "
                                         "try different initial values and/or the option "
                                         "stationary=FALSE");
        }
    }

    arma::mat allProbs(nbObs,nbStates);
    allProbs.ones();

    arma::colvec stepProb(nbObs);
    NumericVector stepArgs(2); // step parameters
    NumericVector step = data["step"];

    arma::colvec angleProb(nbObs);
    NumericVector angleArgs(2); // angle parameters
    NumericVector angle(nbObs);

    if(angleDist!="none") {
        angle = data["angle"];
    }

    arma::rowvec zeromass(nbStates);

    // extract zero-mass parameters from step parameters if necessary
    if(zeroInflation) {
        zeromass = stepPar.row(stepPar.n_rows-1);
        arma::mat stepPar2 = stepPar.submat(0,0,stepPar.n_rows-2,stepPar.n_cols-1);
        stepPar = stepPar2;
    }

    for(int state=0;state<nbStates;state++)
    {
        // compute probabilities of steps

        for(unsigned int i=0;i<stepPar.n_rows;i++)
            stepArgs(i) = stepPar(i,state);

        // if zeroInflation, the probability of zero and non-zero steps must be computed separately
        if(zeroInflation) {
            // remove the NAs from step (impossible to subset a vector with NAs)
            for(int i=0;i<nbObs;i++) {
                if(!std::isfinite(step(i))) {
                    step(i) = -1;
                    stepProb(i) = 1;
                }
            }

            // compute probability of non-zero observations
            stepProb.elem(arma::find(as<arma::vec>(step)>0)) =
                (1-zeromass(state))*funMap[stepDist](step[step>0],stepArgs(0),stepArgs(1));

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
        else {
            stepProb = funMap[stepDist](step,stepArgs(0),stepArgs(1));
        }

        if(angleDist!="none") {
            // compute probabilites of angles
            for(unsigned int i=0;i<anglePar.n_rows;i++)
                angleArgs(i) = anglePar(i,state);

            angleProb = funMap[angleDist](angle,angleArgs(0),angleArgs(1));

            // compute joint probabilities of steps and angles
            allProbs.col(state) = stepProb%angleProb;
        }
        else allProbs.col(state) = stepProb;
    }

    // deal with states known a priori
    double prob = 0;
    if(knownStates(0) != -1) {
        // loop over rows
        for(int i=0 ; i<allProbs.n_rows ; i++) {
            if(knownStates(i)>0) {
                prob = allProbs(i,knownStates(i)-1); // save non-zero probability
                allProbs.row(i).zeros(); // set other probabilities to zero
                allProbs(i,knownStates(i)-1) = prob;
            }
        }
    }

    //======================//
    // 3. Forward algorithm //
    //======================//

    arma::mat Gamma(nbStates,nbStates); // transition probability matrix
    double lscale = 0; // scaled log-likelihood
    int k=0; // animal index
    arma::rowvec alpha(nbStates);

    for(unsigned int i=0;i<allProbs.n_rows;i++) {

        if(nbStates>1)
            Gamma = trMat.slice(i);
        else
            Gamma = 1; // no transition if only one state

        if(k<aInd.size() && i==(unsigned)(aInd(k)-1)) {
            // if 'i' is the 'k'-th element of 'aInd', switch to the next animal
            k++;
            alpha = delta % allProbs.row(i);
        } else {
            alpha = (alpha * Gamma) % allProbs.row(i);
        }

        lscale = lscale + log(sum(alpha));
        alpha = alpha/sum(alpha);
    }

    return -lscale;
}
