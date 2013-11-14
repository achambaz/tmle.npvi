## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
log <- Arguments$getVerbose(-8, timestamp=TRUE);

## set.seed(12345);

sourceDirectory("R");

## Load simulated data set
source("inst/testScripts/simulatedData/01.loadSimulatedDataSet.R")

V <- matrix(runif(3*nrow(obs)), ncol=3)
colnames(V) <- paste("V", 1:3, sep="")
obsV <- cbind(V, obs)  ## 2 covariates !

flavor <- c("learning", "superLearning")[2]

if (TRUE) {
  f <- identity
  ##
  truePsi <- truePsi1
  confInt <- confInt1
  confInt0 <- confInt01
} else {
  f <- function(x) {1*atan(x/1)}
  ##
  truePsi <- truePsi2
  confInt <- confInt2
  confInt0 <- confInt02
}

npvi <- tmle.npvi(obs, f, flavor=flavor, nodes=1)
if (FALSE) {
  npviV <- tmle.npvi(obsV, f, flavor=flavor, nodes=4)
}

if (FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Parameters
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  B <- 1e5+1  ## number of simulated observations for updating estimate of Psi
  tabulate <- TRUE
  family <- "parsimonious"
  
  bound <- 25;
  K <- 10 ## Number of TMLE iterations
  cleverCovTheta <- FALSE  ## use a clever covariate for updating theta ?
  exact <- TRUE
}

history <- getHistory(npvi)
print(round(history, 4))

hp <- history[, "psi"]
hs <- history[, "sic"]
hs[1] <- NA
ics <-  c(-1,1) %*% t(qnorm(0.975)*hs/sqrt(nrow(getObs(npvi))))

## pch <- (20:21)[1:(1+cleverCovTheta)]
pch <- 20
ylim <- range(c(confInt, hp, ics+hp), na.rm=TRUE)

xs <- (1:length(hs))-1
plot(xs, hp, ylim=ylim, pch=pch, xlab="Iteration", ylab=expression(psi[n]))
dummy <- sapply(seq(along=xs), function(x) lines(c(xs[x],xs[x]), hp[x]+ics[, x]))
  
abline(h=confInt, col=4)
abline(h=confInt0, col=2)
  
