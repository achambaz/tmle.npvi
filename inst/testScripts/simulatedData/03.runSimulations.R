## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
log <- Arguments$getVerbose(-8, timestamp=TRUE);
set.seed(12345);

sourceDirectory("R");

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 0. Parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

flavor <- c("learning", "superLearning")[2]
log && cat(log, "Flavor: ", flavor)

useTrueGMu <- FALSE
log && cat(log, "True g and mu: ", useTrueGMu)

B <- 1e5  ## number of simulated observations for updating estimate of Psi
tabulate <- TRUE
family <- "parsimonious"

bound <- 25;
K <- 5 ## Number of TMLE iterations
cleverCovTheta <- TRUE  ## use a clever covariate for updating theta ?
exact <- TRUE

## Learning functions
if (flavor=="superLearning") {
  library(SuperLearner)
  sourceDirectory("inst/testScripts/superLearning")
  source("inst/testScripts/learning/learnCondExpX2givenW.R")  ## no SL equivalent implemented ?
  source("inst/testScripts/learning/learnCondExpXYgivenW.R")  ## no SL equivalent implemented ?

} else {
  sourceDirectory("inst/testScripts/learning");
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 1. Getting true psi (pasted from "01.loadSimulatedDataSet.R")
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lambda0 <- function(W) {
  -W;
}

source("inst/testScripts/realData/01.loadRealDataSet.R")

## Defining class representatives (currently hardcoded)
idxs <- c(NA, 1, 2);  ## row numbers of (up to) 3 class representatives
O <- obs[idxs, ];
O <- as.matrix(O)

## Parameters for the simulation
p <- c(0, 1/2, 1/2);
omega <- c(0, 3, 3); ## initially: c(0, 1, 1)

sCn <- 1.5*(O[3,2]-O[2,2]);
sExp <- (O[3,3]-O[2,3])/4;
## covCnExp <- sqrt(sCn*sExp-1);
covCnExp <- 1;
S <- matrix(c(sCn, covCnExp, covCnExp, sExp), 2 ,2);

sigma2 <- 1

sim0 <- getSample(1e5, O, lambda0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, verbose=log);
str(sim0)
truePsi <- sim0$psi;

verbose && cat(verbose, "'True' psi:");

## True Psi and confidence intervals
cat("\nTrue psi is: ", truePsi, "\n")
confInt0 <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(sim0$obs))
confInt <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(obs))

cat("Optimal confidence interval is: ", confInt, "\n")

t0 <- Sys.time()
profiling <- FALSE
if (profiling) {
  Rprof(tmpfile<-"Rprofiling.txt")
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setting up the bootstrap procedure
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N <- 1e3
nObs <- 2e2

ts <- format(Sys.time(), "%Y-%m-%d-%H-%M-%S")
tag <- sprintf("nObs=%s,N=%d,%s", nObs, N, ts)
if (useTrueGMu) {
  gMuTag <-"trueGMu"
  trueGMu <- list(g=sim0$g, muAux=sim0$muAux)
} else {
  gMuTag <- "estimatedGMu"
  trueGMu <- NULL
}

path <- "simulations"
path <- file.path(path, tolower(flavor), gMuTag, tag)
path <- Arguments$getWritablePath(path)

file1 <- file.path(path, "sim0.RData")

save(sim0, cleverCovTheta, exact, file=file1)

nn <- 0
while (nn<N) {

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Simulation
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dat <- getSample(nObs, O, lambda0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, verbose=log);  
  ## str(dat)
  obs <- dat$obs;

  theX0 <- O[2,2];
  obsC <- obs;
  obsC[, "X"] <- obsC[, "X"] - theX0;
  
  obs <- obsC;
  theX0 <- 0;
  rm(obsC);

  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Declaration
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  c <- 1e-3
  gmin <- 5e-2
  gmax <- 1-gmin
  npvi <- NPVI(obs=obs, family=family, tabulate=tabulate, 
               gmin=gmin, gmax=gmax,
               mumin=min(obs[, "X"])+c, mumax=max(obs[, "X"])-c,
               thetamin=min(obs[, "Y"])+c, thetamax=max(obs[, "Y"])-c)
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Initialization
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bound <- 25
  res <- try(init(npvi, flavor=flavor, 
               learnG=learnG, learnMuAux=learnMuAux, learnTheta=learnTheta,
               bound=bound, B=B, 
               light=TRUE,
               trueGMu=trueGMu, 
               verbose=log))
  if (inherits(res, "try-error")) {
    next;
  }
  
  ## print(getPsi(npvi));
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## A few iterations of TMLE
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for (kk in 1:K) {
    res <- try(update(npvi, flavor=flavor, learnDevG=learnDevG, learnDevMu=learnDevMu, learnDevTheta=learnDevTheta,
                      bound=bound, B=B, cleverCovTheta=cleverCovTheta, exact=exact, useTrueGMu=useTrueGMu, verbose=log));
    if (inherits(res, "try-error")) {
      break;
    }
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Saving results
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  if (!inherits(res, "try-error")) {
    nn <- nn+1
    history <- getHistory(npvi)
    ## print(round(history, 4))
    filename <- paste("history", nn, ".RData", sep="")
    pathname <- file.path(path, filename)
    save(history, file=pathname)
  }
}  
  
if (profiling) {
  Rprof()
  print(summaryRprof(tmpfile)$by.total[1:10,])
}
t1 <- Sys.time()

cat("\n")
print(t1-t0)

