
#' Simulation tool
#'
#' Simulates movement data from an HMM.
#'
#' @param nbAnimals Number of observed individuals to simulate.
#' @param nbStates Number of behavioural states to simulate.
#' @param stepDist Name of the distribution of the step lengths (as a character string).
#' Supported distributions are: gamma, weibull, lnorm, exp. Default: gamma.
#' @param angleDist Name of the distribution of the turning angles (as a character string).
#' Supported distributions are: vm, wrpcauchy. Set to \code{"none"} if the angle distribution should
#' not be estimated. Default: vm.
#' @param stepPar Parameters of the step length distribution.
#' @param anglePar Parameters of the turning angle distribution.
#' @param beta Matrix of regression parameters for the transition probabilities (more information
#' in "Details").
#' @param covs Covariate values to include in the model, as a dataframe. Default: \code{NULL}.
#' Covariates can also be simulated according to a standard normal distribution, by setting
#' \code{covs} to \code{NULL}, and specifying \code{nbCovs>0}.
#' @param nbCovs Number of covariates to simulate (0 by default). Does not need to be specified of
#' \code{covs} is specified.
#' @param zeroInflation \code{TRUE} if the step length distribution is inflated in zero.
#' Default: \code{FALSE}. If \code{TRUE}, values for the zero-mass parameters should be
#' included in \code{stepPar}.
#' @param obsPerAnimal Either the number of the number of observations per animal (if single value),
#' or the bounds of the number of observations per animal (if vector of two values). In the latter case,
#' the numbers of obervations generated for each animal are uniformously picked from this interval.
#' Default: \code{c(500,1500)}.
#' @param model A moveHMM object. This option can be used to simulate from a fitted model.  Default: NULL.
#' Note that, if this argument is specified, most other arguments will be ignored -- except for nbAnimals,
#' obsPerAnimal, covs (if covariate values different from those in the data should be specified),
#' and states.
#' @param states \code{TRUE} if the simulated states should be returned, \code{FALSE} otherwise (default).
#'
#' @return An object moveData, i.e. a dataframe of:
#' \item{ID}{The ID(s) of the observed animal(s)}
#' \item{step}{The step lengths}
#' \item{angle}{The turning angles (if any)}
#' \item{x}{Either easting or longitude}
#' \item{y}{Either northing or latitude}
#' \item{...}{Covariates (if any)}
#'
#' @details \itemize{
#' \item The matrix \code{beta} of regression coefficients for the transition probabilities has
#' one row for the intercept, plus one row for each covariate, and one column for
#' each non-diagonal element of the transition probability matrix. For example, in a 3-state
#' HMM with 2 covariates, the matrix \code{beta} has three rows (intercept + two covariates)
#' and six columns (six non-diagonal elements in the 3x3 transition probability matrix - filled in
#' row-wise).
#' In a covariate-free model (default), \code{beta} has one row, for the intercept.
#'
#' \item If the length of covariate values passed (either through 'covs', or 'model') is not the same
#' as the number of observations suggested by 'nbAnimals' and 'obsPerAnimal', then the series of
#' covariates is either shortened (removing last values - if too long) or extended (starting
#' over from the first values - if too short).
#' }
#'
#' @examples
#' # 1. Pass a fitted model to simulate from
#' # (m is a moveHMM object - as returned by fitHMM - automatically loaded with the package)
#' # We keep the default nbAnimals=1.
#' m <- example$m
#' obsPerAnimal=c(50,100)
#' data <- simData(model=m,obsPerAnimal=obsPerAnimal)
#'
#' # 2. Pass the parameters of the model to simulate from
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' anglePar <- c(pi,0,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                anglePar=anglePar,nbCovs=2,zeroInflation=TRUE,obsPerAnimal=obsPerAnimal)
#'
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(pi,0,0.5,0.7) # mean1, mean2, k1, k2
#' stepDist <- "weibull"
#' angleDist <- "wrpcauchy"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                anglePar=anglePar,obsPerAnimal=obsPerAnimal)
#'
#' # step length only and zero-inflation
#' stepPar <- c(1,10,1,5,0.2,0.3) # mean1, mean2, sd1, sd2, z1, z2
#' stepDist <- "gamma"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist="none",stepPar=stepPar,
#'                nbCovs=2,zeroInflation=TRUE,obsPerAnimal=obsPerAnimal)
#'
#' # include covariates
#' # (note that it is useless to specify "nbCovs", which is determined
#' # by the number of columns of "cov")
#' cov <- data.frame(temp=rnorm(500,20,5))
#' stepPar <- c(1,10,1,5) # mean1, mean2, sd1, sd2
#' anglePar <- c(pi,0,0.5,2) # mean1, mean2, k1, k2
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' data <- simData(nbAnimals=5,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar,
#'                 anglePar=anglePar,covs=cov)
#'
#' @export
#' @importFrom stats rnorm runif

simData <- function(nbAnimals=1,nbStates=2,stepDist=c("gamma","weibull","lnorm","exp"),
                    angleDist=c("vm","wrpcauchy","none"),stepPar=NULL,anglePar=NULL,
                    beta=NULL,covs=NULL,nbCovs=0,zeroInflation=FALSE,obsPerAnimal=c(500,1500),
                    model=NULL,states=FALSE)
{
    ##############################
    ## Check if !is.null(model) ##
    ##############################
    if(!is.null(model)) {
        # extract simulation parameters from model
        nbStates <- ncol(model$mle$stepPar)
        stepDist <- model$conditions$stepDist
        angleDist <- model$conditions$angleDist
        stepPar <- c(t(model$mle$stepPar))
        anglePar <- c(t(model$mle$anglePar))
        beta <- model$mle$beta

        if(is.null(covs)) {
            covs <- model.matrix(model$conditions$formula, model$data)

            if(ncol(covs)>1) {
                # remove intercept column, which is not expected in 'covs'
                names <- colnames(covs)
                covs <- data.frame(covs[,-1]) # data.frame structure is lost when only one column
                colnames(covs) <- names[-1]
            } else
                covs <- NULL
        }
        # else, allow user to enter new values for covariates

        zeroInflation <- model$conditions$zeroInflation

    } else {
        if(is.null(stepPar))
            stop("'stepPar' needs to be specified")

        stepDist <- match.arg(stepDist)
        angleDist <- match.arg(angleDist)
    }

    #####################
    ## Check arguments ##
    #####################
    stepFun <- paste("r",stepDist,sep="")
    angleFun <- paste("r",angleDist,sep="")

    if(nbAnimals<1)
        stop("nbAnimals should be at least 1.")
    if(nbStates<1)
        stop("nbStates should be at least 1.")

    p <- parDef(stepDist,angleDist,nbStates,TRUE,zeroInflation)

    if(length(stepPar)!=p$parSize[1]*nbStates | length(anglePar)!=p$parSize[2]*nbStates) {
        error <- "Wrong number of parameters: there should be"
        error <- paste(error,p$parSize[1]*nbStates,"step parameters and")
        error <- paste(error,p$parSize[2]*nbStates,"angle parameters")
        stop(error)
    }

    if(zeroInflation) {
        stepBounds <- p$bounds[1:((p$parSize[1]-1)*nbStates),]
        sp <- stepPar[1:(length(stepPar)-nbStates)]
        zm <- stepPar[(length(stepPar)-nbStates+1):length(stepPar)]
        if(length(which(zm<0 | zm>1))>0)
            stop("The zero-mass should be in [0,1].")
    } else {
        stepBounds <- p$bounds[1:(p$parSize[1]*nbStates),]
        sp <- stepPar
    }

    if(length(which(sp<=stepBounds[,1] | sp>=stepBounds[,2]))>0)
        stop(paste("Check the step parameters bounds (the parameters should be",
                   "strictly between the bounds of their parameter space)."))

    if(angleDist!="none") {
        # We can't really write distribution-agnostic code here, because the bounds
        # defined in parDef are not the actual bounds of the parameter space.
        m <- anglePar[1:nbStates] # angle mean
        k <- anglePar[(nbStates+1):length(anglePar)] # angle concentration
        if(length(which(m<=(-pi) | m>pi))>0)
            stop("Check the angle parameters bounds. The angle mean should be in (-pi,pi].")
        if(length(which(k<=0))>0)
            stop("Check the angle parameters bounds. The concentration should be strictly positive.")
        if(angleDist=="wrpcauchy" & length(which(k>=1))>0)
            stop("Check the angle parameters bounds. The concentration should be in (0,1).")
    }

    if(length(which(obsPerAnimal<1))>0)
        stop("obsPerAnimal should have positive values.")

    if(!is.null(covs) & nbCovs>0) {
        if(ncol(covs)!=nbCovs)
            warning("covs and nbCovs argument conflicting - nbCovs was set to ncol(covs)")
    }

    if(!is.null(covs)) {
        if(!is.data.frame(covs))
            stop("'covs' should be a data.frame")
    }

    if(!is.null(covs)) {
        nbCovs <- ncol(covs)

        # account for missing values of the covariates
        if(length(which(is.na(covs)))>0)
            warning(paste("There are",length(which(is.na(covs))),
                          "missing covariate values.",
                          "Each will be replaced by the closest available value."))
        for(i in 1:nbCovs) {
            if(length(which(is.na(covs[,i])))>0) { # if covariate i has missing values
                if(is.na(covs[1,i])) { # if the first value of the covariate is missing
                    k <- 1
                    while(is.na(covs[k,i])) k <- k+1
                    for(j in k:2) covs[j-1,i] <- covs[j,i]
                }
                for(j in 2:nrow(trackData))
                    if(is.na(covs[j,i])) covs[j,i] <- covs[j-1,i]
            }
        }
    }

    if(length(obsPerAnimal)==1)
        obsPerAnimal <- rep(obsPerAnimal,2)
    else if(length(obsPerAnimal)!=2)
        stop("obsPerAnimal should be of length 1 or 2.")

    #######################################
    ## Prepare parameters for simulation ##
    #######################################
    # define number of observations for each animal
    allNbObs <- rep(NA,nbAnimals)
    for(zoo in 1:nbAnimals) {
        if(obsPerAnimal[1]!=obsPerAnimal[2])
            allNbObs[zoo] <- sample(obsPerAnimal[1]:obsPerAnimal[2],size=1)
        else
            allNbObs[zoo] <- obsPerAnimal[1]
    }

    # extend covs if not enough covariate values
    if(!is.null(covs)) {
        covnames <- colnames(covs)
        while(sum(allNbObs)>nrow(covs))
            covs <- rbind(covs,covs)
        # shrink covs if too many covariate values
        covs <- data.frame(covs[1:sum(allNbObs),])
        colnames(covs) <- covnames
    }

    # generate regression parameters for transition probabilities
    if(is.null(beta))
        beta <- matrix(rnorm(nbStates*(nbStates-1)*(nbCovs+1)),nrow=nbCovs+1)
    else if(nrow(beta)!=nbCovs+1 | ncol(beta)!=nbStates*(nbStates-1)) {
        if(nbStates>1)
            stop(paste("beta should have ",nbCovs+1," rows and ",nbStates*(nbStates-1)," columns.",sep=""))
        else
            stop("beta should be NULL")
    }

    # initial state distribution
    delta <- rep(1,nbStates)/nbStates

    # format parameters
    wpar <- n2w(c(stepPar,anglePar),p$bounds,beta,delta,nbStates,estAngleMean=(angleDist!="none"))
    par <- w2n(wpar,p$bounds,p$parSize,nbStates,nbCovs,estAngleMean=(angleDist!="none"),stationary=FALSE)

    if(zeroInflation) {
        zeroMass <- par$stepPar[nrow(par$stepPar),]
        stepPar <- par$stepPar[-(nrow(par$stepPar)),,drop=FALSE]
    }
    else {
        zeroMass <- rep(0,nbStates)
        stepPar <- par$stepPar
    }
    anglePar <- par$anglePar # i.e. NULL if angleDist=="none"

    trackData <- NULL
    allCovs <- NULL
    allStates <- NULL

    # build the data frame to be returned
    data <- data.frame(ID=character(),
                       step=numeric(),
                       angle=numeric(),
                       x=numeric(),
                       y=numeric())

    ###########################
    ## Loop over the animals ##
    ###########################
    for (zoo in 1:nbAnimals) {

        # number of observations for animal zoo
        nbObs <- allNbObs[zoo]

        ###############################
        ## Simulate covariate values ##
        ###############################
        if(nbCovs>0) {
            if(is.null(covs)) {
                subCovs <- data.frame(cov1=rnorm(nbObs))
                if(nbCovs>1) {
                    for(j in 2:nbCovs) {
                        c <- data.frame(rnorm(nbObs))
                        colnames(c) <- paste("cov",j,sep="")
                        subCovs <- cbind(subCovs,c)
                    }
                }
            } else {
                # select covariate values which concern the current animal
                if(zoo<2)
                    ind1 <- 1
                else
                    ind1 <- sum(allNbObs[1:(zoo-1)])+1
                ind2 <- sum(allNbObs[1:zoo])
                subCovs <- data.frame(covs[ind1:ind2,])
                if(!is.null(covs))
                    colnames(subCovs) <- colnames(covs) # keep covariates names from input
            }
            allCovs <- rbind(allCovs,subCovs)
        }

        ###############################
        ## Simulate state sequence Z ##
        ###############################
        if(nbStates>1) {
            Z <- rep(NA,nbObs)
            Z[1] <- sample(1:nbStates,size=1,prob=delta)
            for (k in 2:nbObs) {
                gamma <- diag(nbStates)

                g <- beta[1,]
                if(nbCovs==1) g <- g + beta[2,]*subCovs[k,1]
                if(nbCovs>1) {
                    for(j in 1:nbCovs)
                        g <- g + beta[j+1,]*subCovs[k,j]
                }

                gamma[!gamma] <- exp(g)
                gamma <- t(gamma)
                gamma <- gamma/apply(gamma,1,sum)
                Z[k] <- sample(1:nbStates,size=1,prob=gamma[Z[k-1],])
            }
            allStates <- c(allStates,Z)
        } else
            Z <- rep(1,nbObs)

        X <- matrix(nbObs,nrow=nbObs,ncol=2)
        X[1,] <- c(0,0) # initial position of animal

        phi <- 0
        s <- rep(NA,nbObs)
        a <- rep(NA,nbObs)

        ############################
        ## Simulate movement path ##
        ############################
        for (k in 1:(nbObs-1)){
            # prepare lists of arguments for step and angle distributions
            stepArgs <- list(1) ; angleArgs <- list(1) # first argument = 1 (one random draw)
            for(j in 1:nrow(stepPar))
                stepArgs[[j+1]] <- stepPar[j,Z[k]]

            if(angleDist!="none") {
                for(j in 1:nrow(anglePar))
                    angleArgs[[j+1]] <- anglePar[j,Z[k]]
            }

            if(stepDist=="gamma") {
                shape <- stepArgs[[2]]^2/stepArgs[[3]]^2
                scale <- stepArgs[[3]]^2/stepArgs[[2]]
                stepArgs[[2]] <- shape
                stepArgs[[3]] <- 1/scale # rgamma expects rate=1/scale
            }

            if(runif(1)>zeroMass[Z[k]])
                s[k] <- do.call(stepFun,stepArgs)
            else
                s[k] <- 0

            if(angleDist!="none" & s[k]>0) {
                a[k] <- do.call(angleFun,angleArgs)
                if(a[k] >  pi) a[k] <- a[k]-2*pi
                if(a[k] < -pi) a[k] <- a[k]+2*pi
                phi <- phi + a[k]
            }

            m <- s[k]*c(Re(exp(1i*phi)),Im(exp(1i*phi)))
            X[k+1,] <- X[k,] + m
        }

        a[1] <- NA # the first angle value is arbitrary

        # angle = NA if step = 0
        step0 <- which(s==0)
        a[c(step0,step0+1)] <- NA

        d <- data.frame(ID=rep(zoo,nbObs),step=s,angle=a,x=X[,1],y=X[,2])
        data <- rbind(data,d)
    }

    # if covs provided as argument
    if(!is.null(covs) & is.null(allCovs))
        allCovs <- covs

    if(nbCovs>0)
        data <- cbind(data,allCovs)

    # include states sequence in the data
    if(states)
        data <- cbind(data,states=allStates)
    return(moveData(data))
}
