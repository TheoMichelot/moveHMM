
context("deltaMethod")

test_that("Output has the right format",{
  m <- ex$mod
  c <- deltaMethod(m)

  expect_equal(length(c),2)
  expect_equal(dim(c$inf),dim(m$mle$anglePar))
  expect_equal(dim(c$sup),dim(m$mle$anglePar))
})
