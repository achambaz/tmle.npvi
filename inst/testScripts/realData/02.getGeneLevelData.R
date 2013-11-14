library("R.utils")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

## Working on OV data (and the GBM gene EGFR for simulation)

path <- "extdata"
path <- Arguments$getReadablePath(path)
dsName <- "TCGA,OV,expCnMeth"
ts <- "2011-09-23"

for (chr in 1:22) {
  
fileName <- sprintf("%s,chr%s,%s.xdr", dsName, chr, ts)
pathname <- file.path(path, fileName)
## Note for RR: Created by /Users/pierre/Documents/Projects/TCGA/ExpCnMethMirnaData/20100824,dataSets.R

dat <- loadObject(pathname)
str(dat)

geneNames <- dimnames(dat)[["geneSymbol"]];

## finding neutral CN state...
cn <- dat["copyNumber",,]
str(cn)

## moving away from log-ratios
cn <- 2*2^cn
library(aroma.light)

## look across genes for where to cut
vals <- apply(cn, 1, FUN=function(y) {
##  fpv <- findPeaksAndValleys(y, from=min(y), to=4)
  fpv <- findPeaksAndValleys(y)
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

leftVals <- vals[1,][is.finite(vals[1,])]
midPeaks <- vals[2,][is.finite(vals[2,])]
rightVals <- vals[3,][is.finite(vals[3,])]

leftThr <- median(leftVals)
midThr <- median(midPeaks) ## very close to 2.  keep 2.
midThr <- 2
rightThr <- median(rightVals)

## symmetrize !
if (rightThr-midThr < midThr-leftThr) {
  leftThr <- midThr - (rightThr-midThr)
} else {
  rightThr <- midThr + (midThr-leftThr)
}
figPath <- "fig"
figPath <- Arguments$getWritablePath(figPath)
fileName <- sprintf("%s,chr%s,%s.png", dsName, chr, ts)
pathname <- file.path(figPath, fileName)

png(pathname, width=600, height=800)
plotDensity(t(cn[1:50,]), from=0, to=5, main=paste("Chr", chr))
abline(v=c(leftThr, midThr, rightThr))
devDone()
    
centerGeneLevelData <- function(gld) {
  obs <- t(gld)
##  obs <- t(dat[, idx, ])
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

## Save after centering
obsList <- apply(dat, 2, centerGeneLevelData)

path <- "extdata"
path <- Arguments$getWritablePath(path)
fileName <- sprintf("%s,chr%s,%s,centered.xdr", dsName, chr, ts)
pathname <- file.path(path, fileName)

saveObject(obsList, file=pathname)

}
