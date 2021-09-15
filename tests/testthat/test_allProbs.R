
context("allProbs")

test_that("The output has the right format",{
    step <- rgamma(100, 5, 0.25)
    angle <- rvm(100, 0, 2)
    data <- list(step = step, angle = angle)
    stepPar <- matrix(c(8, 20, 5, 10), ncol = 2, byrow = TRUE)
    anglePar <- matrix(c(0, 0, 1.5, 0.7), ncol = 2, byrow = TRUE)
    nbStates <- 2
    p <- allProbs(data = data, nbStates = nbStates, stepDist = "gamma",
                  angleDist = "vm", stepPar = stepPar, anglePar = anglePar)

    expect_equal(nrow(p), length(step))
    expect_equal(ncol(p), nbStates)
})

test_that("It works without turning angles",{
    step <- rgamma(100, 5, 0.25)
    data <- list(step = step)
    stepPar <- matrix(c(8, 20, 5, 10), ncol = 2, byrow = TRUE)
    nbStates <- 2

    expect_error(allProbs(data = data, nbStates = nbStates, stepDist = "gamma",
                          angleDist = "none", stepPar = stepPar), NA)
})

test_that("Zero-inflation works",{
    step <- rgamma(100, 5, 0.25)
    angle <- rvm(100, 0, 2)
    data <- list(step = step, angle = angle)
    stepPar <- matrix(c(8, 20, 5, 10, 0.2, 0.3), ncol = 2, byrow = TRUE)
    anglePar <- matrix(c(0, 0, 1.5, 0.7), ncol = 2, byrow = TRUE)
    nbStates <- 2

    expect_error(allProbs(data = data, nbStates = nbStates, stepDist = "gamma",
                          angleDist = "vm", stepPar = stepPar, anglePar = anglePar,
                          zeroInflation = TRUE), NA)
})
