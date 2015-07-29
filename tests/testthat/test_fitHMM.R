
test_that("Exceptions are thrown",{
  data <- example$data

  nbStates <- 2
  nbCovs <- 2
  stepDist <- "gamma"
  angleDist <- "vm"
  zeroInflation <- FALSE
  angleMean <- c(pi,0)

  # estimation
  mu0 <- c(20,70)
  sigma0 <- c(10,30)
  kappa0 <- c(1,1)
  angleMean <- c(pi,0)
  stepPar0 <- c(mu0,sigma0)
  anglePar0 <- kappa0
  formula <- ~cov1+cos(cov2)
  nbCovs <- length(attr(terms(formula), "term.labels"))

  beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                  nrow=nbCovs+1,byrow=TRUE)
  delta0 <- rep(1,nbStates)/nbStates

  expect_that(fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation), not(throws_error()))

  expect_that(fitHMM(0,data,stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data.frame(),stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data,c(),anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,"unif",angleDist,
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,stepDist,"norm",
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data,-stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation), throws_error())

  expect_that(fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     c(angleMean,0),zeroInflation), throws_error())

})

test_that("The outputs has the right class",{
  data <- example$data

  nbStates <- 2
  nbCovs <- 2
  stepDist <- "gamma"
  angleDist <- "vm"
  zeroInflation <- FALSE
  angleMean <- c(pi,0)

  # estimation
  mu0 <- c(20,70)
  sigma0 <- c(10,30)
  kappa0 <- c(1,1)
  angleMean <- c(pi,0)
  stepPar0 <- c(mu0,sigma0)
  anglePar0 <- kappa0
  formula <- ~cov1+cos(cov2)
  nbCovs <- length(attr(terms(formula), "term.labels"))

  beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                  nrow=nbCovs+1,byrow=TRUE)
  delta0 <- rep(1,nbStates)/nbStates

  mod <- fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,stepDist,angleDist,
                     angleMean,zeroInflation)

  expect_equal(length(which(class(mod)=="moveHMM")),1)
})
