
context("CI")

test_that("Output has the right format",{
  m <- ex$m
  c <- CI(m)

  expect_equal(length(c),2)
  expect_equal(dim(c$lower$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$upper$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$lower$anglePar),dim(m$mle$anglePar))
  expect_equal(dim(c$upper$anglePar),dim(m$mle$anglePar))
  expect_equal(dim(c$lower$beta),dim(m$mle$beta))
  expect_equal(dim(c$upper$beta),dim(m$mle$beta))
})
