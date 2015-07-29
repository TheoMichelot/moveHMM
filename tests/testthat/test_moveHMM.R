
context("moveHMM")

test_that("Exceptions are thrown",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"

  mod <- list(data=data,states=states,mle=mle,stepDist=stepDist,angleDist=angleDist)
  expect_that(moveHMM(mod),not(throws_error()))

  mod <- list(data=data,states=states,mle=mle,stepDist=stepDist)
  expect_that(moveData(data),throws_error())
})

test_that("The output has the right class attribute",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"
  mod <- list(data=data,states=states,mle=mle,stepDist=stepDist,angleDist=angleDist)
  mod <- moveHMM(mod)

  expect_equal(length(which(class(mod)=="moveHMM")),1)
})
