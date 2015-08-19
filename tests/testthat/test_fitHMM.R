
context("fitHMM")

test_that("Exceptions are thrown",{
  data <- ex$data
  simPar <- ex$simPar
  par0 <- ex$par0

  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), not(throws_error()))

  # if nbStates<1
  expect_that(fitHMM(0,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if data empty
  expect_that(fitHMM(simPar$nbStates,data.frame(),par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if stepPar empty
  expect_that(fitHMM(simPar$nbStates,data,c(),par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if stepDist not in list
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,"unif",simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if angleDist not in list
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,"norm",simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if stepPar not within bounds
  expect_that(fitHMM(simPar$nbStates,data,-par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if angleMean too long
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,c(simPar$angleMean,0),
                     simPar$zeroInflation), throws_error())

  # if wrong number of initial parameters
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0[-1],par0$anglePar0,par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0[-1],par0$beta0,par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if beta0 has the wrong dimensions
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0[-1,],par0$delta0,
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

  # if delta0 has the wrong length
  expect_that(fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0[-1],
                     par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                     simPar$zeroInflation), throws_error())

})

test_that("The output has the right class",{
  data <- ex$data
  simPar <- ex$simPar
  par0 <- ex$par0

  mod <- fitHMM(simPar$nbStates,data,par0$stepPar0,par0$anglePar0,par0$beta0,par0$delta0,
                par0$formula,simPar$stepDist,simPar$angleDist,simPar$angleMean,
                simPar$zeroInflation)

  expect_equal(length(which(class(mod)=="moveHMM")),1)
})

test_that("Step length only + zero-inflation works",{
  nbAnimals <- 2
  nbStates <- 2
  nbCovs <- 2
  mu <- c(10,60)
  sigma <- c(10,40)
  zeromass <- c(0.4,0)
  stepPar <- c(mu,sigma,zeromass)
  anglePar <- NULL
  stepDist <- "gamma"
  angleDist <- "none"
  zeroInflation <- TRUE
  nbAnim <- c(50,100)

  data <- simData(nbAnimals=nbAnimals,nbStates=nbStates,stepDist=stepDist,angleDist=angleDist,
                  stepPar=stepPar,anglePar=anglePar,nbCovs=nbCovs,zeroInflation=zeroInflation,
                  obsPerAnimal=nbAnim)

  mu0 <- c(20,50)
  sigma0 <- c(15,30)
  zeromass0 <- c(0.2,0.05)
  stepPar0 <- c(mu0,sigma0,zeromass0)
  anglePar0 <- NULL
  angleMean <- NULL
  formula <- ~cov1+cov2
  beta0 <- NULL
  delta0 <- NULL

  expect_that(fitHMM(nbStates,data,stepPar0,anglePar0,beta0,delta0,formula,
              stepDist,angleDist,angleMean,zeroInflation,verbose=0),
              not(throws_error()))
})
