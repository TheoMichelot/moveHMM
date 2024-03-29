
context("turnAngle")

test_that("expected values",  {
    a <- c(0, 0)
    b <- c(0, 1)
    c <- c(1, 1)
    d <- c(0, 2)
    expect_equal(turnAngle(a, b, c, FALSE), -pi/2)
    expect_equal(turnAngle(a, b, a, FALSE), pi)
    expect_equal(turnAngle(a, b, d, FALSE), 0)
})
