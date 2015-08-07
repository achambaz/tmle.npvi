##
## testing H20
##

library(tmle.npvi)

set.seed(12345)

log <- Arguments$getVerbose(-20, timestamp=TRUE)

## generating a data set

O <- cbind(W=c(0.05218652, 0.01113460),
           X=c(2.722713, 9.362432),
           Y=c(-0.4569579, 1.2470822))
O <- rbind(NA, O)
lambda0 <- function(W) {-W}
p <- c(0, 1/2, 1/2)
omega <- c(0, 3, 3)
S <- matrix(c(10, 1, 1, 0.5), 2 ,2)

sim <- getSample(2e2, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
obs <- sim$obs
V <- matrix(runif(3*nrow(obs)), ncol=3)
colnames(V) <- paste("V", 1:3, sep="")
obs <- cbind(V, obs)
X0 <- O[2,2]
obsC <- obs
obsC[, "X"] <- obsC[, "X"] - X0
obs <- obsC

## 

if (FALSE) {
  library(SuperLearner)
  ## works only when 'nodes' equals 1 without the previous command...
  nodes <- 2
  npvi.SL <- tmle.npvi(obs, f=identity,
                       flavor="superLearning",
                       nodes=nodes, cvControl=list(V=2),
                       B=5e4, nMax=10, 
                       verbose=log)
}

if (TRUE) {
  ## problem: need to do that first:
  library(h2o)
  library(h2oEnsemble)
  ## debug(h2oEnsemble:::.make_Z)
  nodes <- 3
  npvi.EL <- tmle.npvi:::tmle.npvi.(obs, f=identity,
                                    flavor="h2oEnsembleLearning",
                                    nodes=nodes, cvControl=list(V=2),
                                    B=5e4, nMax=10, 
                                    verbose=log)
}
