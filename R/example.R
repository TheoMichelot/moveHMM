
#' Example data simulation
#'
#' Simulates the file data/example.RData, used in other functions examples.

example <- function()
{
  # simulate data
  nbAnimals <- 1
  nbStates <- 2
  nbCovs <- 2
  mu<-c(15,50)
  sigma<-c(10,20)
  angleMean <- c(pi,0)
  kappa <- c(0.7,1.5)
  stepPar <- c(mu,sigma)
  anglePar <- c(angleMean,kappa)
  stepDist <- "gamma"
  angleDist <- "vm"
  zeroInflation <- FALSE
  data <- simData(nbAnimals,nbStates,stepDist,angleDist,stepPar,anglePar,nbCovs,zeroInflation)

  # estimation
  mu0 <- c(20,70)
  sigma0 <- c(10,30)
  kappa0 <- c(1,1)
  stepPar0 <- c(mu0,sigma0)
  anglePar0 <- kappa0
  formula <- ~cov1+cos(cov2)
  nbCovs <- length(attr(terms(formula), "term.labels"))

  beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                  nrow=nbCovs+1,byrow=TRUE)
  delta0 <- rep(1,nbStates)/nbStates

  mod <- fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,
                "gamma","vm",angleMean,zeroInflation)

  example <- list(data=data,mod=mod)
  save(example,file="data/example.RData")
}
