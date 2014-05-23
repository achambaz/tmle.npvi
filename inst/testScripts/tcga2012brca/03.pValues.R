getPValue <- function( ## Calculates p-value from an object of type 'history'
                      history,
### The \code{history} of a TMLE procedure.
                      nobs
### An \code{integer}, the associated number of observations.
                      ){
  ##seealso<< tmle.npvi, getHistory, as.character.NPVI
  
  y <- history[nrow(history), ]
  psi <- y["psi"]
  phi <- y["phi"]
  se <- y["sicAlt"]/sqrt(nobs)
  pval <- 2*pnorm(abs(psi-phi), sd=se, lower.tail=FALSE)
  names(pval) <- "p.value"
  return(pval)
### Returns the p-value of the two-sided test of ``\eqn{Psi(P_0)=Phi(P_0)}''.
}


PVAL <- sapply(tmle, function(ll){ifelse(!is.null(ll$hist),
                                         getPValue(ll[[2]], 463), NA)})
yi <- -log10(PVAL[!is.na(PVAL)])
pos <- sapply(names(yi), function(ll){rev(unlist(strsplit(ll, split=",")))[1]})
pos <- 1e-6*as.integer(pos)
geneNames <- sapply(names(yi), function(ll){unlist(strsplit(ll, split=","))[2]})

dev.new()
thr <- 1.5
ww <- which(yi>thr)

ylim <- c(0, max(yi))
rg <- range(pos)
xlim <- rg*c(.95, 1.05)

##png(pathname, width=width, height=height)
par(cex=2, mar=c(5, 4, 0, 0)+.2)
plot(NA, xlim=xlim, ylim=ylim,
     xlab="Genome position (Mb)", ylab="Test statistic")
abline(h=thr, col=2)
pusr <- par()$usr
unitX <-.01*(pusr[2]-pusr[1])
unitY <-.02*(pusr[4]-pusr[3])
pchs <- rep(20, length=length(pos))
if (length(ww)) {
  pchs[ww] <- 1
  if (require(maptools)) {
    pointLabel(pos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=.4) # col=cols[ww], )
  } else {
    text(pos[ww]+unitX, yi[ww]+unitY, labels=geneNames[ww], cex=0.5) # col=cols[ww])
  }
}
points(pos[-ww], yi[-ww], pch=pchs[-ww], cex=1)



