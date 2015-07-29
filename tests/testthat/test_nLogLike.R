
context("nLogLike")

test_that("Exceptions are thrown",{
  d <- example$data

  nbStates <- 2
  nbCovs <- 2
  stepDist <- "gamma"
  angleDist <- "vm"
  mu0 <- c(20,80)
  sigma0 <- c(20,40)
  kappa0 <- c(1,1)
  beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                  nrow=nbCovs+1,byrow=TRUE)
  delta0 <- c(1,1)/2
  par0 <- c(mu0,sigma0,kappa0)

  bounds <- parDef(stepDist,angleDist,nbStates,FALSE,FALSE)$bounds
  parSize <- parDef(stepDist,angleDist,nbStates,FALSE,FALSE)$parSize

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)
  angleMean <- c(pi,0)

  expect_that(nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,FALSE),
              not(throws_error()))

  expect_that(nLogLike(wpar[-1],nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,FALSE),
              throws_error())

  expect_that(nLogLike(wpar,nbStates,bounds,parSize,data,"unif",angleDist,angleMean,FALSE),
              throws_error())

  expect_that(nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,"norm",angleMean,FALSE),
              throws_error())

  data <- data[-2] # remove data$step
  expect_that(nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,FALSE),
              throws_error())
})

test_that("angleMean=NULL, angleDist=NULL, and zeroInflation=TRUE work",{
  d <- example$data

  nbStates <- 2
  nbCovs <- 2
  stepDist <- "gamma"
  angleDist <- "NULL"
  mu0 <- c(20,80)
  sigma0 <- c(20,40)
  kappa0 <- c(1,1)
  beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                  nrow=nbCovs+1,byrow=TRUE)
  delta0 <- c(1,1)/2
  par0 <- c(mu0,sigma0,kappa0)

  bounds <- parDef(stepDist,angleDist,nbStates,TRUE,TRUE)$bounds
  parSize <- parDef(stepDist,angleDist,nbStates,TRUE,TRUE)$parSize

  wpar <- n2w(par0,bounds,beta0,delta0,nbStates)
  angleMean <- NULL

  expect_that(nLogLike(wpar,nbStates,bounds,parSize,data,stepDist,angleDist,angleMean,TRUE),
              not(throws_error()))
})
