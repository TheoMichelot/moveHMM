
context("moveHMM")

test_that("Exceptions are thrown",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"
  mod <- list(1,2,3)

  m <- list(data=data,states=states,mle=mle,stepDist=stepDist,angleDist=angleDist,mod=mod,
            zeroInflation=FALSE)
  expect_that(moveHMM(m),not(throws_error()))

  m <- list(data=data,states=states,mle=mle,stepDist=stepDist)
  expect_that(moveHMM(m),throws_error())
})

test_that("The output has the right class attribute",{
  data <- data.frame(1,2,3)
  states <- rep(1,10)
  mle <- list(1,2,3)
  stepDist <- "gamma"
  angleDist <- "vm"
  mod <- list(1,2,3)

  m <- list(data=data,states=states,mle=mle,stepDist=stepDist,angleDist=angleDist,mod=mod,
            zeroInflation=FALSE)
  m <- moveHMM(m)

  expect_equal(length(which(class(m)=="moveHMM")),1)
})
