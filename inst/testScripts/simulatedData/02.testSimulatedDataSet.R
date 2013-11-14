## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
log <- Arguments$getVerbose(-8, timestamp=TRUE);

## set.seed(12345);

sourceDirectory("R");

## Load simulated data set
## source("inst/testScripts/simulatedData/01.loadSimulatedDataSet.R")

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

flavor <- c("learning", "superLearning")[2]
B <- 1e5+1  ## number of simulated observations for updating estimate of Psi
useTrueGMu <- FALSE
if (useTrueGMu) {
  trueGMu <- list(g=sim0$g, muAux=sim0$muAux)
} else {
  trueGMu <- NULL
}

tabulate <- TRUE
family <- "parsimonious"

bound <- 25;
K <- 4 ## Number of TMLE iterations
cleverCovTheta <- FALSE  ## use a clever covariate for updating theta ?
exact <- TRUE

## Learning functions
if (flavor=="superLearning") {
  library(SuperLearner)
  sourceDirectory("inst/testScripts/superLearning")
} else {
  sourceDirectory("inst/testScripts/learning");
}

t0 <- Sys.time()
profiling <- FALSE
if (profiling) {
  Rprof(tmpfile<-"Rprofiling.txt")
}


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Declaration
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
npvi <- NPVI(obs=obs, family=family, tabulate=tabulate, 
             gmin=5e-2, gmax=.95,
             mumin=min(obs[, "X"]), mumax=max(obs[, "X"]),
             thetamin=min(obs[, "Y"]), thetamax=max(obs[, "Y"]))

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Initialization
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

init(npvi, flavor=flavor, 
     learnG=learnG, learnMuAux=learnMuAux, learnTheta=learnTheta,
     bound=bound, B=B,
     light=TRUE,
     trueGMu=trueGMu, 
     verbose=log);
print(getPsi(npvi));



for (kk in 1:K) {
  update(npvi, flavor=flavor, learnDevG=learnDevG, learnDevMu=learnDevMu, learnDevTheta=learnDevTheta,
         bound=bound, B=B, cleverCovTheta=cleverCovTheta, exact=exact, trueGMu=trueGMu, verbose=log);
  print(kk)
}

history <- getHistory(npvi)
print(round(history, 4))

hp <- history[, "psi"]
hs <- history[, "sic"]
hs[1] <- NA
ics <-  c(-1,1) %*% t(qnorm(0.975)*hs/sqrt(nrow(getObs(npvi))))

pch <- (20:21)[1:(1+cleverCovTheta)]
ylim <- range(c(confInt, hp, ics+hp), na.rm=TRUE)

xs <- (1:length(hs))-1
plot(xs, hp, ylim=ylim, pch=pch, xlab="Iteration", ylab=expression(psi[n]))
dummy <- sapply(seq(along=xs), function(x) lines(c(xs[x],xs[x]), hp[x]+ics[, x]))

abline(h=confInt, col=4)
abline(h=confInt0, col=2)

## psi <- sapply(res, getPsi)
## psiPn <- sapply(res, getPsiPn)
## mic1 <- sapply(res[-1], function(xx){mean(getEfficientInfluenceCurve(xx)[, 1])})
## mic2 <- sapply(res[-1], function(xx){mean(getEfficientInfluenceCurve(xx)[, 2])})
## mic <- sapply(res[-1], function(xx){mean(getEfficientInfluenceCurve(xx)[, 3])})
## eps <- sapply(res[-1], getEpsilon)
## logLikIncr <- sapply(res[-1], getLogLikIncr)

## impact <- evalImpact(res, obs)


if (profiling) {
  Rprof()
  print(summaryRprof(tmpfile)$by.total[1:10,])
}
t1 <- Sys.time()

cat("\n")
print(t1-t0)
