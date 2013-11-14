library("R.utils")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

## GBM not found (probably never actually created..)
## Working on OV instead (and the GBM gene EGFR for simulation)

pathname <- "~/Documents/Projects/TargetedMLE/data/TCGA,OV,expCnMeth,2010-08-24.xdr"
## Note for RR: Created by /Users/pierre/Documents/Projects/TCGA/ExpCnMethMirnaData/20100824,dataSets.R

dat <- loadObject(pathname)
str(dat)

geneNames <- dimnames(dat)[["geneSymbol"]];

## finding neutral CN state...
cn <- dat["copyNumber",,]
str(cn)

## looks like data are already centered at CN=0 !
plot(density(apply(cn, 1, median)))

## moving away from log-ratios
cn <- 2*2^cn
library(aroma.light)
plotDensity(t(cn[1:10,]), from=0, to=4)

## where to cut ?
fpv <- findPeaksAndValleys(cn[1,], from=0, to=4)
fpv <- findPeaksAndValleys(cn[2,], from=0, to=4)
fpv <- subset(fpv, density>1e-2)
## closest peak to CN=2
fp <- subset(fpv, type=="peak")
neutral <- which.min(abs(fp$x-2))

## closest valleys
left <- max(subset(fpv, type=="valley" & x<2)$x)
right <- min(subset(fpv, type=="valley" & x>2)$x)

## look across genes
vals <- apply(cn, 1, FUN=function(y) {
  fpv <- findPeaksAndValleys(y, from=0, to=4)
  fpv <- subset(fpv, density>1e-2)
  ## closest peak to CN=2
  fp <- subset(fpv, type=="peak")
  neutral <- which.min(abs(fp$x-2))
  mid <- fp$x[neutral]
  ## valleys closest to CN=2
  left <- max(subset(fpv, type=="valley" & x<2)$x)
  right <- min(subset(fpv, type=="valley" & x>2)$x)
  c(left, mid, right)
})

sum(!is.finite(vals))
plotDensity(t(vals))  ## sounds good

leftVals <- vals[1,][is.finite(vals[1,])]
midPeaks <- vals[2,][is.finite(vals[2,])]
rightVals <- vals[3,][is.finite(vals[3,])]

boxplot(leftVals, midPeaks, rightVals)  ## pretty tight: OK

leftThr <- median(leftVals)
midThr <- median(midPeaks) ## very close to 2.  keep 2.
midThr <- 2
rightThr <- median(rightVals)

getObs <- function(idx) {
  obs <- t(dat[, idx, ])
  colnames(obs) <- c("W", "X", "Y")

  cc <- which(complete.cases(obs))
  obsC <- obs[cc, ]

  theX0 <- midThr
  X <- 2*2^obsC[, "X"]
  idx <- which(leftThr<=X & X<=rightThr)
  obsC[idx, "X"] <- theX0;
  obsC[idx, "X"] <- obsC[idx, "X"] - theX0;
  theX0 <- 0;

  obsC;
}

o <- getObs(2)
pairs(o)
## pairs(o[,c("Y", "X", "W")])

## a gene
geneList <- c("TP53", "NF1", "BRCA1", "BRCA2", "RB1", "CDK12", "CCNE1", "FOXM1")
geneList <- paste("^", geneList, "$", sep="")
geneList <- c(geneList, "NOTCH")
ggList <- unlist(lapply(geneList, grep, geneNames))
names(ggList) <- geneNames[ggList]

ggList

for (ii in seq(along=ggList)) {
  o <- getObs(ggList[ii])
  pairs(o, main=names(ggList[ii]))
  grDevices::devAskNewPage(ask = TRUE)
}

