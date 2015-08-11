
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
