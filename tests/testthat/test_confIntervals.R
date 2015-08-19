
context("confIntervals")

test_that("Output has the right format",{
  m <- ex$mod
  c <- confIntervals(m)

  expect_equal(length(c),2)
  expect_equal(dim(c$inf$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$sup$stepPar),dim(m$mle$stepPar))
  expect_equal(dim(c$inf$beta),dim(m$mle$beta))
  expect_equal(dim(c$sup$beta),dim(m$mle$beta))
})

test_that("inf<estimate<sup",{
  m <- ex$mod
  c <- confIntervals(m)

  # check that inferior bound is lower than the estimate
  expect_equal(length(which(c$inf$stepPar<=m$mle$stepPar)),length(m$mle$stepPar))
  expect_equal(length(which(c$inf$beta<=m$mle$beta)),length(m$mle$beta))

  # check that superior bound is higher than the estimate
  expect_equal(length(which(c$sup$stepPar>=m$mle$stepPar)),length(m$mle$stepPar))
  expect_equal(length(which(c$sup$beta>=m$mle$beta)),length(m$mle$beta))
})
