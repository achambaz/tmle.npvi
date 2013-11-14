## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
library(SuperLearner)

log <- Arguments$getVerbose(-8, timestamp=TRUE);

set.seed(12345);

sourceDirectory("R");

path <- "extdata"
path <- Arguments$getReadablePath(path)
fileName <- "TCGA,OV,expCnMeth,2010-08-24,centered,2011-09-20.xdr"
pathname <- file.path(path, fileName)
obsList <- loadObject(pathname)
geneNames <- names(obsList)

## Learning functions
sourceDirectory("inst/testScripts/superLearning");
source("inst/testScripts/learning/learnCondExpX2givenW.R")  ## no SL equivalent implemented ?
source("inst/testScripts/learning/learnCondExpXYgivenW.R")  ## no SL equivalent implemented ?

flavor <- "superLearning";

B <- 1e5
  
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA") ##antoine## , "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");

learnTheta <- c("SL.glm.theta", "SL.polymars", SL.library);
## learnTheta0 <- c("SL.glm.theta", "SL.polymars", SL.library);  ## not tested 
learnG <- c("SL.glm.g", SL.library);
## learnMuAux <- c("SL.glm.muAux", SL.library)
learnMuAux <- c(SL.library);

## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA", "SL.glmnet");
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm") ##antoine##, "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");

learnDevG <- c(SL.library)
learnDevMu <- c(SL.library)
learnDevTheta <- c(SL.library)

useTrueGMu <- FALSE
tabulate <- TRUE
family <- "parsimonious"

t0 <- Sys.time()


bound <- 0.1;  ## more than 0.1 would not really make sense

nGenes <- 6261
nbrIt <- 5; ## number of TMLE iterations
resMat <- matrix(NA, nGenes, 1+nbrIt);
colnames(resMat) <- sprintf("psi%s", as.character(0:nbrIt));

epsMat <- matrix(NA, nGenes, nbrIt);
colnames(epsMat) <- sprintf("epsilon%s", 1:nbrIt);

varMat <- matrix(NA, nGenes, 1+nbrIt);
colnames(varMat) <- sprintf("var%s", 0:nbrIt);

## for (gg in 1:nGenes) {
for (gg in 6) {
  geneName <- geneNames[gg]
  log && enter(log, sprintf("Gene %s/%s: %s", gg, nGenes, geneName));
  obs <- obsList[[gg]]
  str(obs)

##  obs <- as.data.frame(obs)

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
  bound <- 25;
  log && enter(log, "Initialization")
  tt <- try(init(npvi, flavor=flavor,
                 learnG=learnG, learnMuAux=learnMuAux, learnTheta=learnTheta,
                 bound=bound, B=B,
                 light=TRUE,
                 useTrueGMu=useTrueGMu, 
                 verbose=more(log, 0)))
  log && exit(log)
  if (class(tt)!="try-error") {
    resMat[gg, "psi0"] <- getPsi(npvi)
    effIC <- getEfficientInfluenceCurve(npvi)
    varIC <- var(effIC[,"eic"]);
    ## confInt <- getPsi(npvi)+c(-1, 1)*qnorm(.975)*sqrt(varIC/nrow(obs))      
    varMat[gg, "var0"] <- varIC
  } else {
    next
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cleverCovTheta <- TRUE
  exact <- TRUE
  
  for (kk in 1:nbrIt) {
    print(kk)
    tt <- try(update(npvi, flavor=flavor, learnDevG=learnDevG, learnDevMu=learnDevMu,
                     learnDevTheta=learnDevTheta,
                     bound=bound, B=B, cleverCovTheta=cleverCovTheta,
                     exact=exact, useTrueGMu=useTrueGMu, verbose=log))
    if (class(tt)!="try-error") {
      print(getPsi(npvi))
      resMat[gg, sprintf("psi%s", kk)] <- getPsi(npvi)
      epsMat[gg, sprintf("epsilon%s", kk)] <- getEpsilon(npvi)
      ## need to be able to estimate  a confidence interval as well at each iteration !
      effIC <- getEfficientInfluenceCurve(npvi)
      varIC <- var(effIC[,"eic"]);
      varMat[gg, sprintf("var%s", kk)] <- varIC
      ## confInt <- getPsi(npvi)+c(-1, 1)*qnorm(.975)*sqrt(varIC/nrow(obs))      
    } else {
      break
    }
  }
  print(resMat[gg,])
}
psi0 <- resMat[, "psi0"];
psi1 <- resMat[, "psi1"];

## boxplot(resMat);
## plot(psi0-psi, psi1-psi);
## plot(psi0, psi1)
## abline(a=0, b=1, col=2)

