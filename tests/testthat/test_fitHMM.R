
context("fitHMM")

test_that("Exceptions are thrown", {
    data <- example$data
    simPar <- example$simPar
    par0 <- example$par0

    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean),  NA)

    # if nbStates<1
    expect_error(fitHMM(data = data, nbStates = 0, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if data empty
    expect_error(fitHMM(data = data.frame(), nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if stepPar empty
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = c(), anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if stepDist not in list
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = "unif",
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if angleDist not in list
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = "norm", angleMean = simPar$angleMean))

    # if stepPar not within bounds
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = -par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if angleMean too long
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = c(simPar$angleMean, 0)))

    # if wrong number of initial parameters
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0[-1], anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0[-1],
                       beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if beta0 has the wrong dimensions
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0[1, ], delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

    # if delta0 has the wrong length
    expect_error(fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                       beta0 = par0$beta0, delta0 = par0$delta0[-1], formula = par0$formula, stepDist = simPar$stepDist,
                       angleDist = simPar$angleDist, angleMean = simPar$angleMean))

})

test_that("The output has the right class", {
    data <- example$data
    simPar <- example$simPar
    par0 <- example$par0

    m <- fitHMM(data = data, nbStates = simPar$nbStates, stepPar0 = par0$stepPar0, anglePar0 = par0$anglePar0,
                beta0 = par0$beta0, delta0 = par0$delta0, formula = par0$formula, stepDist = simPar$stepDist,
                angleDist = simPar$angleDist, angleMean = simPar$angleMean)

    expect_equal(length(which(class(m) == "moveHMM")), 1)
})

test_that("Step length only + zero-inflation works", {
    set.seed(1)
    nbAnimals <- 2
    nbStates <- 2
    nbCovs <- 2
    mu <- c(10, 60)
    sigma <- c(10, 40)
    zeromass <- c(0.4, 0)
    stepPar <- c(mu, sigma, zeromass)
    anglePar <- NULL
    stepDist <- "gamma"
    angleDist <- "none"
    zeroInflation <- TRUE
    nbAnim <- c(50, 100)

    data <- simData(nbAnimals = nbAnimals, nbStates = nbStates, stepDist = stepDist, angleDist = angleDist,
                    stepPar = stepPar, anglePar = anglePar, nbCovs = nbCovs, zeroInflation = zeroInflation,
                    obsPerAnimal = nbAnim)

    mu0 <- c(10, 60)
    sigma0 <- c(10, 40)
    zeromass0 <- c(0.4, 0.01)
    stepPar0 <- c(mu0, sigma0, zeromass0)
    anglePar0 <- NULL
    angleMean <- NULL
    formula <- ~cov1+cov2

    expect_error(fitHMM(data = data, nbStates = nbStates, stepPar0 = stepPar0, anglePar0 = anglePar0, formula = formula,
                       stepDist = stepDist, angleDist = angleDist, angleMean = angleMean, verbose = 0),  NA)
})
