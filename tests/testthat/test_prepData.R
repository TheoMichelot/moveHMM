
context("prepData")

test_that("Exception is thrown", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    y <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2)
    z <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2)

    trackData <- data.frame(x)
    expect_error(prepData(trackData))

    trackData <- data.frame(x, z)
    expect_error(prepData(trackData))

    trackData <- data.frame(x, y)
    expect_error(prepData(trackData),  NA)
})

test_that("The right slots are defined", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    y <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2)
    trackData <- data.frame(x, y)
    data <- prepData(trackData)

    expect_true(!is.null(data$ID))
    expect_true(!is.null(data$x))
    expect_true(!is.null(data$y))
    expect_true(!is.null(data$step))
    expect_true(!is.null(data$angle))
})

test_that("The returned object is of the correct class", {
    x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    y <- c(1, 1, 1, 2, 2, 2, 1, 1, 1, 2)
    trackData <- data.frame(x, y)
    data <- prepData(trackData)

    expect_equal(class(data), c("moveData", "data.frame"))
})
