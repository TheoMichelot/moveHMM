
context("viterbi")

test_that("Exceptions are thrown",{
  data <- example$data
  mod <- example$mod
  nbStates <- 2
  stepDist <- "gamma"
  angleDist <- "vm"
  angleMean <- c(pi,0)

  # reconstruction of states sequence
  expect_that(viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar,
                    mod$mle$anglePar,angleMean),not(throws_error()))

  expect_that(viterbi(data.frame(),nbStates,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar,
                      mod$mle$anglePar,angleMean),throws_error())

  expect_that(viterbi(data,0,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar,
                      mod$mle$anglePar,angleMean),throws_error())

  expect_that(viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,"unif",angleDist,mod$mle$stepPar,
                      mod$mle$anglePar,angleMean),throws_error())

  expect_that(viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,stepDist,"norm",mod$mle$stepPar,
                      mod$mle$anglePar,angleMean),throws_error())

  expect_that(viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar[-1],
                      mod$mle$anglePar,angleMean),throws_error())

  expect_that(viterbi(data,nbStates,mod$mle$beta,mod$mle$delta,stepDist,angleDist,mod$mle$stepPar,
                      mod$mle$anglePar,c(angleMean,1)),throws_error())
})
