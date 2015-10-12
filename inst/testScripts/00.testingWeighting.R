##
## testing weighting
##

library(tmle.npvi)

set.seed(12345)

## generating a data set

O <- cbind(W=c(0.05218652, 0.01113460),
           X=c(2.722713, 9.362432),
           Y=c(-0.4569579, 1.2470822))
O <- rbind(NA, O)
lambda0 <- function(W) {-W}
p <- c(0, 1/2, 1/2)
omega <- c(0, 3, 3)
S <- matrix(c(10, 1, 1, 0.5), 2 ,2)

sampleSize <- 1e4

sim <- getSample(sampleSize, O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
obs <- sim$obs
V <- matrix(runif(3*nrow(obs)), ncol=3)
colnames(V) <- paste("V", 1:3, sep="")
obs <- cbind(V, obs)
X0 <- O[2,2]
obsC <- obs
obsC[, "X"] <- obsC[, "X"] - X0
obs <- obsC

## randomly removing observations with X!=0

isZero <- (obs[,"X"]==0)
probs <- rep(1, nrow(obs))
probs[!isZero] <- 1/2
keep <- rbinom(length(probs), 1, probs)
obs2 <- obs[keep==1, ]

## preparing weights

weights <- rep(1, nrow(obs2))
weights[obs2[, "X"]!=0] <- 2
weights <- weights/sum(weights)

## inferring psi

npvi2 <- tmle.npvi(obs2, f=identity, weights=weights, flavor="learning", B=5e4, nMax=10)

## comparing with estimator based on as many observations as in 'obs2'

sim <- getSample(nrow(obs2), O, lambda0, p=p, omega=omega, sigma2=1, Sigma3=S)
obs <- sim$obs
V <- matrix(runif(3*nrow(obs)), ncol=3)
colnames(V) <- paste("V", 1:3, sep="")
obs <- cbind(V, obs)
X0 <- O[2,2]
obsC <- obs
obsC[, "X"] <- obsC[, "X"] - X0
obs3 <- obsC

npvi3 <- tmle.npvi(obs3, f=identity, flavor="learning", B=5e4, nMax=10)
