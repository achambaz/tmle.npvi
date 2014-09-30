getPValue <- function(# Calculates p-value from an object of type 'history'
                      history,
### The \code{history} of a TMLE procedure.
                      nobs,
### An \code{integer}, the associated number of observations.
                      wrt.phi=TRUE
### A  \code{logical}  equal  to  \code{TRUE}  by default,  which  means  that
### \eqn{psi_n}  is  compared  with  \eqn{phi_n}.  Otherwise,  \eqn{psi_n}  is
### compared with 0.
){
  ##seealso<< tmle.npvi, getHistory, as.character.NPVI
  y <- history[nrow(history), ]
  psi <- y["psi"]
  if (wrt.phi) {
    phi <- y["phi"]
    se <- y["sicAlt"]/sqrt(nobs)
  } else {
    phi <- 0
    se <- y["psi.sd"]/sqrt(nobs)
  }
  pval <- 2*pnorm(abs(psi-phi), sd=se, lower.tail=FALSE)
  names(pval) <- "p.value"
  return(pval)
### Returns the p-value of the two-sided test of ``\eqn{Psi(P_0)=Phi(P_0)}''.
}

what <- "WholeGenome"
path <- Arguments$getReadablePath(system.file("testScripts/tcga2012brca/",
                                              package="tmle.npvi"))
path <- "."

PVAL <- NULL
for (flavor in c("learning", "superLearning")) {
  pathname <- file.path(path, paste(flavor, what, ".xdr", sep=""))
  tmle <- loadObject(pathname)$TMLE
  if (flavor=="learning") {
    tmle.l <- tmle
  } else {
    tmle.sl <- tmle
  }
  pathname <- file.path(path, "cumLimChr.xdr")
  cumLimChr <- loadObject(pathname)
  
  pval <- sapply(tmle, function(ll){getPValue(ll$hist, 463, TRUE)})
  PVAL[[flavor]] <- pval
  rm(tmle)

  yi <- -log10(pval)
  if (any(is.infinite(yi))) {
    yi[is.infinite(yi)] <- max(yi[!is.infinite(yi)])+1 ## arbitrary
  }
  chr <- sapply(names(yi), function(ll){unlist(strsplit(ll, split=","))[1]})
  chr <- sapply(chr, function(ll){unlist(strsplit(ll, split="chr"))[2]})
  chr <- as.integer(chr)
  posRel <- sapply(names(yi), function(ll){unlist(strsplit(ll, split=","))[2]})
  posRel <- as.integer(posRel)*1e-3
  posAbs <- posRel + cumLimChr[chr]
  geneNames <- sapply(names(yi), function(ll){unlist(strsplit(ll, split=","))[3]})
  geneNames <- sapply(geneNames, function(ll){unlist(strsplit(ll, split="\\."))[1]})
  attributes(geneNames) <- NULL

  dev.new()
  thr <- 200
  ww <- which(yi>thr)
  
  ylim <- c(0, 300) ## c(0, max(yi))
  rg <- range(posAbs)
  xlim <- rg*c(.95, 1.05)
  
  ##png(pathname, width=width, height=height)
  par(cex=2, mar=c(5, 4, 2, 0)+.2)
  plot(NA, xlim=xlim, ylim=ylim,
       xlab=paste("Genome position\n-", what, "-"),
       ylab="-log10(pval)",
       main=flavor, axes=FALSE)
  abline(h=thr, col=2)
  pusr <- par()$usr
  unitX <-.01*(pusr[2]-pusr[1])
  unitY <-.02*(pusr[4]-pusr[3])
  pchs <- rep(20, length=length(posAbs))
  if (length(ww)) {
    pchs[ww] <- .5
    if (require(maptools)) {
      points(posAbs[ww], yi[ww], pch=pchs[ww], cex=0.25)
      pointLabel(posAbs[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww], )
    } else {
      text(posAbs[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww])
    }
  }
  points(posAbs[-ww], yi[-ww], pch=pchs[-ww], cex=0.5)
  abline(v=cumLimChr, col="orange")
  xx <- sapply(1:(length(cumLimChr)-1), function(ii) {
    mean(cumLimChr[0:1+ii])
  })
  box()
  axis(2)
  axis(1, xx, 1:length(xx), tcl=NA)
}

cmp <- sapply(1:length(tmle.l), function(ii){c(tmle.l[[ii]]$hist[nrow(tmle.l[[ii]]$hist), "psi"],
                                               tmle.sl[[ii]]$hist[nrow(tmle.sl[[ii]]$hist), "psi"])})


stop()



ww <- apply(sapply(PVAL, is.na), 1, any)
ranks <- matrix(c(order(PVAL[[1]][!ww]),
                  order(PVAL[[2]][!ww])),
                nrow=2, byrow=TRUE)
geneNames <- sapply(names(PVAL[[1]][!ww]), function(ll){unlist(strsplit(ll, split=","))[2]})
attributes(geneNames) <- NULL
rownames(ranks) <- names(PVAL)

library(RankAggreg)
cmprRanks <- RankAggreg(ranks, 10, method="CE", distance="Spearman",
                        N=100, convIn=5, rho=.1, verbose=FALSE)
top <- geneNames[cmprRanks$top.list]
