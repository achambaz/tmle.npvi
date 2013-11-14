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
dsName <- "TCGA,OV,expCnMeth"
chr <- 18
ts <- "2011-09-23"

fileName <- sprintf("%s,chr%s,%s,centered.xdr", dsName, chr, ts)
pathname <- file.path(path, fileName)
obsList <- loadObject(pathname)
geneNames <- names(obsList)

flavor <- "superLearning";
## Learning functions
sourceDirectory("inst/testScripts/superLearning");
## source("inst/testScripts/learning/learnCondExpX2givenW.R")  ## no SL equivalent implemented ?
## source("inst/testScripts/learning/learnCondExpXYgivenW.R")  ## no SL equivalent implemented ?


B <- 1e5
tabulate <- TRUE
family <- "parsimonious"

t0 <- Sys.time()

bound <- 0.1;  ## more than 0.1 would not really make sense

nGenes <- length(obsList)
print(nGenes)

nbrIt <- 5; ## number of TMLE iterations
resMat <- matrix(NA, nGenes, 1+nbrIt);
colnames(resMat) <- sprintf("psi%s", as.character(0:nbrIt));
rownames(resMat) <- geneNames
  
epsMat <- matrix(NA, nGenes, nbrIt);
colnames(epsMat) <- sprintf("epsilon%s", 1:nbrIt);
rownames(epsMat) <- geneNames

varMat <- matrix(NA, nGenes, 1+nbrIt);
colnames(varMat) <- sprintf("var%s", 0:nbrIt);
rownames(varMat) <- geneNames

span <- 1:nGenes
##span <- 1:100 + 400
##span <- 501:nGenes
##span <- 1:65+65
span <- seq(from=94, to=nGenes, by=1)
span <- 58

print(span)
path <- "results"
path <- Arguments$getWritablePath(path)
ts <- format(Sys.Date())
fileName <- sprintf("%s,chr%s,%s,%s.rda", dsName, chr, paste(range(span), collapse="-"), ts)
pathname <- file.path(path, fileName)

for (gg in span) {
  geneName <- geneNames[gg]
  log && enter(log, sprintf("Gene %s/%s: %s", gg, nGenes, geneName));
  obs <- obsList[[gg]]
  str(obs)

##  obs <- as.data.frame(obs)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Declaration
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  c <- 1e-3
  npvi <- NPVI(obs=obs, family=family, tabulate=tabulate, 
               gmin=5e-2, gmax=.95,
               mumin=min(obs[, "X"])+c, mumax=max(obs[, "X"])-c,
               thetamin=min(obs[, "Y"])+c, thetamax=max(obs[, "Y"])-c)
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Initialization
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bound <- 25;
  log && enter(log, "Initialization")
  tt <- try(init(npvi, flavor=flavor,
                 learnG=learnG, learnMuAux=learnMuAux, learnTheta=learnTheta,
                 bound=bound, B=B,
                 light=TRUE,
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
                     exact=exact, verbose=log))
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
  save(resMat, epsMat, varMat, file=pathname)
}
psi0 <- resMat[, "psi0"];
psi1 <- resMat[, "psi1"];

## boxplot(resMat);
## plot(psi0-psi, psi1-psi);
## plot(psi0, psi1)
## abline(a=0, b=1, col=2)


if (FALSE) {
  res <- tail(getHistory(npvi), 1)[1, ]
##  res <- getHistory(npvi)[4, ]

  ## testing "psi_0 = 0"
  res[["psi"]]
  sqrt(nrow(obs))*res[["psi"]]/res[["sic"]]

  ## testing "psi_0 = phi_0"
  res[["phi"]]
  sqrt(nrow(obs))*(res[["psi"]]-res[["phi"]])/res[["sicAlt"]]
}
