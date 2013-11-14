## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Setup
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

library(R.utils);
library(SuperLearner)

log <- Arguments$getVerbose(-8, timestamp=TRUE);

nulls <- c("zero", "correlation coefficient", "partial correlation coefficient", "linear regression coefficient")

thresholds <- c(45, 20, 20, 6)  ## threshold for the statistic to be plotted
names(thresholds) <- nulls

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

nObss <- sapply(obsList, nrow)  ## different numbers of observations per gene
summary.factor(nObss)
## end: load original data

## load results
path <- "realData"
path <- Arguments$getReadablePath(path)

patt <- sprintf("%s,chr%s,.*.rda", dsName, chr)
fileNames <- list.files(path, pattern=patt)

fileName <- fileNames[1] ## yopa
pathname <- file.path(path, fileName)
load(pathname)
str(resMat)

pn <- "extdata/TCGA,OV,expCnMeth,chr18,annotation,2011-09-29.xdr"
if (file.exists(pn)) {
  annDat <- loadObject(pn)
} else if (require(biomaRt)) {
  source("inst/testScripts/realData/20110925,getGenePositions.R")
} else {
  throw("Could not retrieve annotations either from file:" , pn, " or using biomaRt\nPlease install biomaRt")
}

pos <- (annDat[, "start"] + annDat[, "end"])/2

cor.coef <- sapply(obsList, FUN=function(obs) {
  cor(obs[, "X"], obs[, "Y"])
})

cor.p <- sapply(obsList, FUN=function(obs) {
  cor.test(obs[, "X"], obs[, "Y"])$p.value
})

cor.stat <- sapply(obsList, FUN=function(obs) {
  cor.test(obs[, "X"], obs[, "Y"])$statistic
})

parCor.coef <- sapply(obsList, FUN=function(obs) {
  pcor(obs)[3,2]
})

linReg.coef <- sapply(obsList, FUN=function(obs) {
  XX <- obs[, "X"]
  YY <- obs[, "Y"]
  sum(XX*YY)/sum(XX*XX)
})


## p-values: what is H0 ??
zeros <- rep(0, length(geneNames))
names(zeros) <- geneNames
nullValues <- list("zero"=zeros,
                    "correlation coefficient"=cor.coef,
                    "partial correlation coefficient"=parCor.coef,
                    "linear regression coefficient"=linReg.coef
                    )
nullValues <- nullValues[nulls]  ## in case the ordering would be different
str(nullValues)

nObsMat <- matrix(nObss, nrow(resMat), ncol(resMat))
nullMatList <- lapply(nullValues, matrix, nrow(resMat), ncol(resMat))
str(nullMatList)

zMatList <- lapply(nullValues, FUN=function(nullValue) {
  sqrt(nObsMat)*(resMat-nullValue)/sqrt(varMat)
})
str(zMatList)

pMatList <- lapply(zMatList, FUN=function(zMat) 1-pnorm(zMat)) ## get p-values
str(pMatList)

qMatList <- lapply(pMatList, p.adjust, method="BH") ## multiple testing correction
lapply(qMatList, FUN=function(qMat) print(sum(qMat==0, na.rm=TRUE)))

## which k to choose ?

signifList <- lapply(seq(along=nulls), FUN=function(nn) {
  signif<- apply(zMatList[[nn]], 2, FUN=function(x) which(abs(x)>thresholds[[nn]]))
})
str(signifList)

chosenK <- 3  ## 3rd TMLE iteration
## idxs <- 1+c(0, chosenK)
## cols <- c("darkgray", "black")
idxs <- 1+chosenK  ## plot only the chosen one
cols <- 1

figPath <- "png"
figPath <- Arguments$getWritablePath(figPath)

width <- 800
height <- 800

## which should be colored differently ?
w1 <- signifList[[1]][[idxs]]
w2 <- signifList[[4]][[idxs]]
cols <- rep(1, length(geneNames))
wint <- intersect(w1, w2)
cols[wint] <- "blue"

            
for (nn in seq(along=nulls)) {
  null <- nulls[[nn]]
  nullTag <- toCamelCase(null)
  thr <- thresholds[[null]]
  figName <- sprintf("%s,chr%s,%s,|z|>%s", dsName, chr, nullTag, thr)
  filename <- sprintf("%s.png", figName)
  pathname <- file.path(figPath, filename)

  ## if (null %in% c("zero", "linear regression coefficient")) {
  datMat <- zMatList[[null]]
  ylab <- "Test statistic"
  ## } else {
  ##   datMat <- -log10(qMatList[[null]])
  ##   ylab <- "-log10(q-value)"
  ## }

  ylim <- range(datMat, na.rm=TRUE)
  pos <- 1e-6*(annDat[, "start"] + annDat[, "end"])/2
  rg <- range(pos, na.rm=TRUE)
  xlim <- rg*c(.95, 1.05)

  png(pathname, width=width, height=height)
  par(cex=2, mar=c(5, 4, 0, 0)+.2)
  plot(NA, xlim=xlim, ylim=ylim, xlab="Genome position (Mb)", ylab=ylab)
  abline(h=0)
  pusr <- par()$usr
  unitX <-.01*(pusr[2]-pusr[1])
  unitY <-.02*(pusr[4]-pusr[3])
  for (ii in seq(along=idxs)) {  
    yi <- datMat[, idxs[ii]]
    ww <- which(abs(yi)>thr)
    pchs <- rep(20, length=length(geneNames))
    if (length(ww)) {
      pchs[ww] <- 1
      if (require(maptools)) {
        pointLabel(pos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], col=cols[ww], cex=.8)
      } else {
        text(pos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], col=cols[ww], cex=0.5)
      }
    }
    points(pos, yi, col=cols, pch=pchs, cex=1)
  }
  devDone()

  signif <- signifList[[null]]
  ## - plot original data for genes of interest
  for (ww in signif[[chosenK+1]]) {
    geneName <- geneNames[ww]
    figName <- sprintf("%s,chr%s,%s,|z|>%s,%s", dsName, chr, nullTag, thr, geneName)
    filename <- sprintf("%s.png", figName)
    pathname <- file.path(figPath, filename)
    
    png(pathname, width=width, height=height)
    par(cex=2, mar=c(5, 4, 2, 0)+.2)
    dat <- obsList[[ww]]
    colnames(dat) <- c("DNA methylation", "DNA copy number", "gene expression")
    pairs(dat, main=geneName)
    devDone()
  }
}
  
  
