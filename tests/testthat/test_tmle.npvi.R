library("tmle.npvi")

context("Flavors of the function tmle.npvi")

set.seed(12345)
##
## Simulating a data set and computing the true value of the parameter
##

## Parameters for the simulation (case 'f=identity')
O <- cbind(W=c(0.05218652, 0.01113460),
           X=c(2.722713, 9.362432),
           Y=c(-0.4569579, 1.2470822))
O <- rbind(NA, O)
lambda0 <- function(W) {-W}
p <- c(0, 1/2, 1/2)
omega <- c(0, 3, 3)
S <- matrix(c(10, 1, 1, 0.5), 2 ,2)

## Simulating a data set of 200 i.i.d. observations
sim <- getSample(2e2, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
obs <- sim$obs

## Adding (dummy) baseline covariates
V <- matrix(runif(3*nrow(obs)), ncol=3)
colnames(V) <- paste("V", 1:3, sep="")
obs <- cbind(V, obs)

## Caution! MAKING '0' THE REFERENCE VALUE FOR 'X'
obs[, "X"] <- obs[, "X"] - O[2,2]

test_that("Flavor 'learning' gives expected results", {
    expect_equal_to_reference(tmle.npvi(obs, f=identity, flavor="learning", B=5e4, nMax=10),
                              "learning.rds")
})

test_that("Flavor 'superLearning' gives expected results", {
    expect_equal_to_reference(tmle.npvi(obs, f=identity, flavor="superLearning", B=5e4, nMax=10),
                              "superLearning.rds")
})

test_that("Flavor 'h2oEnsembleLearning' gives expected results", {
#    skip("Not working currently--update in h2o/h2oEnsemble?")
    npvi <- tmle.npvi(obs, f=identity, flavor="h2oEnsembleLearning", B=5e4, nMax=10)
#    expect_equal_to_reference(tmle.npvi(obs, f=identity, flavor="h2oEnsembleLearning", B=5e4, nMax=10),
#                              "h2oEnsembleLearning.rds")
})

