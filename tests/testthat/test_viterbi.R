
context("viterbi")

test_that("Output has the right format",{
  m <- ex$m
  expect_equal(length(viterbi(m)),nrow(m$data))
})
