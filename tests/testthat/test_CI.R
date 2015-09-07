
context("CI")

test_that("Output has the right format",{
  m <- ex$mod
  c <- CI(m)

  expect_equal(length(c),2)
  expect_equal(dim(c$inf$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$sup$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$inf$anglePar),dim(m$mle$anglePar))
  expect_equal(dim(c$sup$anglePar),dim(m$mle$anglePar))
  expect_equal(dim(c$inf$beta),dim(m$mle$beta))
  expect_equal(dim(c$sup$beta),dim(m$mle$beta))
})
