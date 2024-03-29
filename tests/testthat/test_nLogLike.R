
context("nLogLike")

test_that("Exceptions are thrown", {
    data <- example$data
    simPar <- example$simPar
    par0 <- example$par0

    estAngleMean <- is.null(simPar$angleMean)
    bounds <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                     estAngleMean, simPar$zeroInflation)$bounds
    parSize <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                      estAngleMean, simPar$zeroInflation)$parSize

    par <- c(par0$stepPar0, par0$anglePar0)
    wpar <- n2w(par, bounds, par0$beta0, par0$delta0, simPar$nbStates, estAngleMean)

    expect_error(nLogLike(wpar, simPar$nbStates, bounds, parSize, data, simPar$stepDist, simPar$angleDist,
                         simPar$angleMean, simPar$zeroInflation), NA)

    # if not enough parameters provided
    expect_error(nLogLike(wpar[-1], simPar$nbStates, bounds, parSize, data, simPar$stepDist, simPar$angleDist,
                         simPar$angleMean, simPar$zeroInflation))

    # if stepDist not in list
    expect_error(nLogLike(wpar, simPar$nbStates, bounds, parSize, data, "unif", simPar$angleDist,
                         simPar$angleMean, simPar$zeroInflation))

    # if angleDist not in list
    expect_error(nLogLike(wpar, simPar$nbStates, bounds, parSize, data, simPar$stepDist, "norm",
                         simPar$angleMean, simPar$zeroInflation))

    data <- data[-2] # remove data$step
    expect_error(nLogLike(wpar, simPar$nbStates, bounds, parSize, data, simPar$stepDist, simPar$angleDist,
                         simPar$angleMean, simPar$zeroInflation))
})

test_that("angleMean = NULL, angleDist = NULL, and zeroInflation = TRUE work", {
    data <- example$data
    simPar <- example$simPar
    par0 <- example$par0

    estAngleMean <- TRUE
    bounds <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                     estAngleMean, TRUE)$bounds
    parSize <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                      estAngleMean, TRUE)$parSize

    par0$stepPar0 <- c(par0$stepPar0, rep(0.2, simPar$nbStates)) # include zero mass parameters
    par0$anglePar0 <- c(rep(0, simPar$nbStates), par0$anglePar0) # include angle mean parameters
    par <- c(par0$stepPar0, par0$anglePar0)
    wpar <- n2w(par, bounds, par0$beta0, par0$delta0, simPar$nbStates, estAngleMean)

    expect_error(nLogLike(wpar, simPar$nbStates, bounds, parSize, data, simPar$stepDist, "none",
                         NULL, TRUE), NA)
})

test_that("logAlpha, logBeta, and nLogLike are consistent", {
    m <- example$m
    data <- m$data
    simPar <- example$simPar
    nbAnimals <- simPar$nbAnimals
    par <- c(t(c(m$mle$stepPar)), c(t(m$mle$anglePar))[3:4])

    estAngleMean <- m$conditions$estAngleMean
    bounds <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                     estAngleMean, simPar$zeroInflation)$bounds
    parSize <- parDef(simPar$stepDist, simPar$angleDist, simPar$nbStates,
                      estAngleMean, simPar$zeroInflation)$parSize

    nll <- nLogLike(m$mod$estimate, simPar$nbStates, bounds, parSize, data, simPar$stepDist, simPar$angleDist,
                    simPar$angleMean, simPar$zeroInflation)

    la<-logAlpha(m)
    lb<-logBeta(m)
    ll<-0
    for(i in 1:nbAnimals){
        aInd <- max(which(data$ID == i))
        c <- max(la[aInd, ]+lb[aInd, ]) # cancels below ; prevents numerical errors
        ll <- ll + c + log(sum(exp(la[aInd, ]+lb[aInd, ]-c)))
    }
    expect_equal(nll, -ll)

    # random time step from each individual
    for(i in 1:nbAnimals){
        data <- example$m$data
        aInd <- which(data$ID == i)

        data <- data[aInd, ]
        nll <- nLogLike(m$mod$estimate, simPar$nbStates, bounds, parSize, data, simPar$stepDist, simPar$angleDist,
                      simPar$angleMean, simPar$zeroInflation)

        samp <- sample(aInd, 1)
        c <- max(la[samp, ]+lb[samp, ]) # cancels below ; prevents numerical errors
        ll <- c + log(sum(exp(la[samp, ]+lb[samp, ]-c)))
        expect_equal(nll, -ll)
    }
})
