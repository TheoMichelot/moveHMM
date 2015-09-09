
context("viterbi")

test_that("Exceptions are thrown",{
  data <- ex$data
  mod <- ex$m
  simPar <- ex$simPar

  expect_that(viterbi(data,simPar$nbStates,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      simPar$angleDist,mod$mle$stepPar,mod$mle$anglePar,simPar$angleMean),
              not(throws_error()))

  # if data is empty
  expect_that(viterbi(data.frame(),simPar$nbStates,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      simPar$angleDist,mod$mle$stepPar,mod$mle$anglePar,simPar$angleMean),
              throws_error())

  # if nbStates < 1
  expect_that(viterbi(data,0,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      simPar$angleDist,mod$mle$stepPar,mod$mle$anglePar,simPar$angleMean),
              throws_error())

  # if stepDist not in list
  expect_that(viterbi(data,simPar$nbStates,mod$mle$beta,mod$mle$delta,"unif",
                      simPar$angleDist,mod$mle$stepPar,mod$mle$anglePar,simPar$angleMean),
              throws_error())

  # if angleDist not in list
  expect_that(viterbi(data,simPar$nbStates,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      "norm",mod$mle$stepPar,mod$mle$anglePar,simPar$angleMean),
              throws_error())

  # if not enough parameters
  expect_that(viterbi(data,simPar$nbStates,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      simPar$angleDist,mod$mle$stepPar[-1],mod$mle$anglePar,simPar$angleMean),
              throws_error())

  # if angleMean has the wrong size
  expect_that(viterbi(data,simPar$nbStates,mod$mle$beta,mod$mle$delta,simPar$stepDist,
                      simPar$angleDist,mod$mle$stepPar,mod$mle$anglePar,c(simPar$angleMean,1)),
              throws_error())
})
