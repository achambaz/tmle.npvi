## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
log <- Arguments$getVerbose(-8, timestamp=TRUE);

## set.seed(12345678);

sourceDirectory("R");


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Getting true psi
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

source("inst/testScripts/realData/01.loadRealDataSet.R")

theta0 <- function(W) {
  -W;
}

## Defining class representatives (currently hardcoded)

idxs <- c(NA, 1, 2);  ## row numbers of (up to) 3 class representatives
O <- obs[idxs, ];


## Parameters for the simulation
p <- c(0, 1/2, 1/2);
omega <- c(0, 3, 3); ## initially: c(0, 1, 1)


sCn <- 1.5*(O[3,2]-O[2,2]);
sExp <- (O[3,3]-O[2,3])/4;
## covCnExp <- sqrt(sCn*sExp-1);
covCnExp <- 1;
S <- matrix(c(sCn, covCnExp, covCnExp, sExp), 2 ,2);

sigma2 <- 1

sim0 <- getSample(1e5, O, theta0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, verbose=log);
str(sim0)
truePsi <- sim0$psi;

verbose && cat(verbose, "'True' psi:");

cat("\nTrue psi is: ", truePsi, "\n")
confInt0 <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(sim0$obs))
confInt <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(obs))

cat("Optimal confidence interval is: ", confInt, "\n")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Learning functions
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sourceDirectory("inst/testScripts/learning");

flavor <- "learning";

B <- 1e4

useTrueGMu <- FALSE
## family <- "parsimonious"
family <- "gaussian"

t0 <- Sys.time()
profiling <- FALSE
if (profiling) {
  Rprof(tmpfile<-"Rprofiling.txt")
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setting up the bootstrap procedure
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N <- 1e3-380

## folder <- file.path("results", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"))
## dir.create(folder)
folder <- "results/2011-03-17-10-11-13/"
file1 <- file.path(folder, "sim0.RData")

cleverCovTheta <- TRUE
exact <- TRUE

save(sim0, cleverCovTheta, exact, file=file1)

K <- 5
bound <- 25;

nn <- 0
while (nn<N) {

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Simulation
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  dat <- getSample(1e3, O, theta0, p=p, omega=omega, sigma2=sigma2, Sigma3=S, verbose=log);  
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
  npvi <- NPVI(
               gmin=5e-2, gmax=.95,
               mumin=min(obs[, "X"]), mumax=max(obs[, "X"]),
               thetamin=min(obs[, "Y"]), thetamax=max(obs[, "Y"]))
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Initialization
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  init(npvi, obs, flavor=flavor, 
       learnG=learnG, learnMuAux=learnMuAux, learnTheta=learnTheta,
       bound=bound, B=B,
       light=TRUE,
       useTrueGMu=useTrueGMu, 
       verbose=log);
  
  ## print(getPsi(npvi));
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## A few iterations of TMLE
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  ## Remark: o if (a) 'cleverCovTheta' is set to FALSE
  ##              (b) 'learnTheta' outputs the true theta
  ##              (c) 'learnDevTheta' outputs 0
  ##           then TMLE knows and uses the true theta!
  ##         o if 'useTrueGMu' is set to TRUE
  ##           then TMLE knows and uses the true mu and g!
  
  
  for (kk in 1:K) {
    res <- try(update(npvi, obs, flavor=flavor, learnDevG=learnDevG, learnDevMu=learnDevMu, learnDevTheta=learnDevTheta,
                      bound=bound, B=B, leverCovTheta=cleverCovTheta, exact=exact, useTrueGMu=useTrueGMu, verbose=log));
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
    
    save(history, file=file.path(folder, paste("history", nn+380, ".RData", sep="")))
  }
}  
  
if (profiling) {
  Rprof()
  print(summaryRprof(tmpfile)$by.total[1:10,])
}
t1 <- Sys.time()

cat("\n")
print(t1-t0)
  

