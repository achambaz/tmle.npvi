## If given the true g, theta, mu, sigma2, do we recover the true Psi ?

library(R.utils);
log <- Arguments$getVerbose(-8, timestamp=TRUE);

set.seed(2);


sourceDirectory("R");
source("inst/testScripts/simulatedData/01.loadSimulatedDataSet.R")

cat("\nTrue psi is: ", truePsi, "\n")
confInt0 <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(sim0$obs))
confInt <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(dat$varIC/nrow(obs))
cat("Optimal confidence interval is: ", confInt, "\n")

## centering sim0$obs...
obs <- sim0$obs;

theX0 <- O[2,2];
obsC <- obs;
obsC[, "X"] <- obsC[, "X"] - theX0;

obs <- obsC;
theX0 <- 0;
rm(obsC);


## Taking "true" learning functions to get a new sample
trueG <- sim0$g
trueMu <- sim0$mu
trueTheta <- sim0$theta
trueTheta0 <- sim0$theta0
trueSigma2 <- mean(obs[, "X"]^2)

obsB <- simulateData(B, obs[, "W"], obs[, "X"], trueG, trueMu, trueSigma2, family="gaussian")
res0 <- estimatePsi(theta=trueTheta, theta0=trueTheta0, obs=obsB, obsX=obsB[, "X"], verbose=verbose);
psi0 <- res0$mean
print(truePsi)
print(psi0)
## Is this expected/compatible with the theory ?

## sanity check
res0 <- estimatePsi(theta=trueTheta, theta0=trueTheta0, obs=obs, obsX=obs[, "X"], verbose=verbose);
psi0 <- res0$mean
stopifnot(psi0==truePsi)
