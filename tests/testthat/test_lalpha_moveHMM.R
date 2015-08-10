
context("lalpha.moveHMM")

test_that("Exceptions are thrown",{
  m <- example$mod
  expect_that(lalpha(m),not(throws_error()))

  expect_that(lalpha(1),throws_error())
})

test_that("Output has the right format",{
  m <- example$mod
  la <- lalpha(m)

  expect_equal(nrow(la),nrow(m$data))
  expect_equal(ncol(la),ncol(m$mle$stepPar))
})
