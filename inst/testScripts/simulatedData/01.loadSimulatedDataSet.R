## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library("R.utils");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

sourceDirectory("R")  ## to be removed when package is stable

lambda0 <- function(W) {
  -W;
}

source("inst/testScripts/realData/01.loadRealDataSet.R")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Defining class representatives (currently hardcoded)
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
idxs <- c(NA, 1, 2);  ## row numbers of (up to) 3 class representatives
O <- obs[idxs, ];
O <- as.matrix(O)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Parameters for the simulation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
p <- c(0, 1/2, 1/2);
omega <- c(0, 3, 3); ## initially: c(0, 1, 1)

  
sCn <- 1.5*(O[3,2]-O[2,2]);
sExp <- (O[3,3]-O[2,3])/4;
## covCnExp <- sqrt(sCn*sExp-1);
covCnExp <- 1;
S <- matrix(c(sCn, covCnExp, covCnExp, sExp), 2 ,2);

sigma2 <- 1

N <- 2e2

dat <- getSample(N, O, lambda0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, verbose=log);  
str(dat)
obs <- dat$obs;

theX0 <- O[2,2];
obsC <- obs;
obsC[, "X"] <- obsC[, "X"] - theX0;

obs <- obsC;
theX0 <- 0;
rm(obsC);
    
## True Psi and confidence intervals: identity
f1 <- identity
sim01 <- getSample(1e3, O, lambda0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, f=f1, verbose=log);
truePsi1 <- sim01$psi;

cat("\nCase f=identity:\n")
cat("\ttrue psi is: ", signif(truePsi1, 3), "\n")
confInt01 <- truePsi1+c(-1, 1)*qnorm(.975)*sqrt(sim01$varIC/nrow(sim01$obs))
confInt1 <- truePsi1+c(-1, 1)*qnorm(.975)*sqrt(sim01$varIC/nrow(obs))
cat("\toptimal confidence interval is: ", signif(confInt1, 3), "\n")

## True Psi and confidence intervals: atan
f2 <- function(x) {1*atan(x/1)}
sim02 <- getSample(1e3, O, lambda0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, f=f2, verbose=log);
truePsi2 <- sim02$psi;

cat("\nCase f=atan:\n")
cat("\ttrue psi is: ", signif(truePsi2, 3), "\n")
confInt02 <- truePsi2+c(-1, 1)*qnorm(.975)*sqrt(sim02$varIC/nrow(sim02$obs))
confInt2 <- truePsi2+c(-1, 1)*qnorm(.975)*sqrt(sim02$varIC/nrow(obs))
cat("\toptimal confidence interval is: ", signif(confInt2, 3), "\n")

############################################################################
# HISTORY:
# 2010-07-18
# o Created.
############################################################################


