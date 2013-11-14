## - - - - - - - - - - - - - - - - - - - - - - - - 
## general script for the analysis of simulations
## - - - - - - - - - - - - - - - - - - - - - - - - 

anderson.darling.test <- function(x, var)
{
  ## Anderson-Darling test of normality for estimated mean, known variance
  ##
  ## see http://en.wikipedia.org/wiki/Anderson%E2%80%93Darling_test
  if (missing(var)) {
    throw("Missing 'var' in call to anderson.darling.test\n")
  }
  x <- (x-mean(x))/sqrt(var)
  x <- sort(x)
  stat <- -length(x) - mean(seq(1, 2*length(x)-1, by=2)*(log(pnorm(x)) + log(1-pnorm(rev(x)))))
  return(stat)
}

## .libPaths(".")
library(R.utils)
log <- Arguments$getVerbose(-8, timestamp=TRUE);

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 0. Parameters
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## choose only one if desired...
flavors <- c("learning", "superLearning")[1:2][2]
useTrueGMus <- c(TRUE, FALSE)[1:2][2]
doGraphs <- c(TRUE, FALSE)[2]

for (flavor in flavors) {
  for (useTrueGMu in useTrueGMus) {
    
    log && cat(log, "Flavor: ", flavor)
    log && cat(log, "True g and mu: ", useTrueGMu)

    N <- 1e3
    nObs <- 2e2

    tag <- sprintf("nObs=%s,N=%d,", nObs, N)
    if (useTrueGMu) {
      gMuTag <-"trueGMu"
    } else {
      gMuTag <- "estimatedGMu"
    }

    path <- "simulations"
    path <- file.path(path, tolower(flavor), gMuTag)
    path <- Arguments$getReadablePath(path)

    lf <- list.files(path, pattern=tag)
    if (length(lf)==0) {
      stop("No available results for this set of parameters")
    } else {
      ## Taking the first folder
      path <-  file.path(path, lf[1])
    }
    file1 <- file.path(path, "sim0.RData")
    load(file1)

    filename <- "historyList.RData"
    pnl <- file.path(path, filename)
    if (file.exists(pnl)) {
      load(pnl) ## retrieves 'NPVI', a list of 'history'-like matrices
    } else {
      ## creates 'NPVI' from files
      patt <- "history[0-9]+.RData"
      lf <- list.files(path, pattern=patt)
      NPVI <- vector("list", length(lf))
      for (ii in 1:length(lf)) {
        filename <- paste("history", ii, ".RData", sep="")
        pni <- file.path(path, filename)
        load(pni)
        NPVI[[ii]] <- history
      }
      if (length(NPVI)==N) {
        save(NPVI, file=pnl)
      }
    }

##    NPVI <- NPVI[1:100]
    
    
    ## retrieve estimates of psi at each *end* of iteration (useful when cleverCovTheta=TRUE, because then two estimates of psi are calculated at each iteration)
    history <- NPVI[[1]]
    allSteps <- as.numeric(gsub("step([0-9]+)", "\\1", row.names(history)))
    steps <- sapply(unique(allSteps), FUN=function(x) max(which(allSteps==x)))
    ## steps <- seq(from=1, to=nrow(history), by=2)  ## a bit simpler, but less safe...

    ## only look at 3 TMLE iterations
    steps <- steps[1:4]

    figPath <- "png"
    figPath <- Arguments$getWritablePath(figPath)
    figTag <- sprintf("%s,%s,N=%d",  tolower(flavor), gMuTag, N)

    ## empirical densities
    truePsi <- sim0$psi
    confInt0 <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nrow(sim0$obs))
    xlim <- truePsi+c(-2, 2)*qnorm(.975)*sqrt(sim0$varIC/nObs)
    psi <- sapply(NPVI, function(mat){mat[, "psi"]})
    sic <- sapply(NPVI, function(mat){mat[, "sic"]})

    l1 <- apply(abs(psi-truePsi), 1, mean)
    l2 <- sqrt(apply((psi-truePsi)^2, 1, mean))

    densList <- apply(psi, 1, density, from=xlim[1], to=xlim[2])
    densMax <- sapply(densList, function(x) max(x$y))
    ylim <- c(0, max(densMax))
    ylim <- c(0, 30)
    xlim <- range(densList[[1]]$x)

    ttl <- expression("Density of "*psi[n]^k)
    ttl <- ""

    width <- 800
    height <- 800

    figName <-"densities"
    filename <- sprintf("%s,%s.png", figName, figTag)
    pathname <- file.path(figPath, filename)

    if (doGraphs) {
      png(pathname, width=width, height=height)
    } else {
      x11()
    }
    par(cex=3, lwd=3, mar=c(2, 2, 0, 0)+.2)
    plot(NA, xlab="", ylab="", xlim=xlim, ylim=ylim, main=ttl, t='n')
    usr <- par("usr")
    
    confInt <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nObs)
    
    rect(confInt[1], usr[3], confInt[2], usr[4], col="#CCCCCC", border=NA)
    rect(confInt0[1], usr[3], confInt0[2], usr[4], col="#AAAAAA", border=NA)
    box()

    for (ii in seq(along=steps)) {
      step <- steps[ii]
      lines(densList[[step]], col=ii)
    }
    legend("topright", paste("k=", seq(along=steps)-1), col=seq(along=steps), lwd=2)
    if (doGraphs) {
      dev.off()
    }

    figName <-"densities0vs1"
    filename <- sprintf("%s,%s.png", figName, figTag)
    pathname <- file.path(figPath, filename)

    if (doGraphs) {
      png(pathname, width=width, height=height)
    } else {
      x11()
    }
    par(cex=3, lwd=1.5, mar=c(2, 2, 2, 0)+.2)
    d0 <- density(psi[steps[1], ], from=xlim[1], to=xlim[2])
    d1 <- density(psi[steps[2], ], from=xlim[1], to=xlim[2])
    plot(psi[steps[1],], psi[steps[2],], main=expression(psi[n]^1*" vs "*psi[n]^0), pch=20, cex=0.2, xlim=xlim, ylim=xlim)
    usr <- par("usr")
    rect(confInt0[1], usr[3], confInt0[2], usr[4], col="gray", border=1)
    usr <- par("usr")
    rect( usr[1], confInt0[1], usr[2], confInt0[2], col="gray", border=2)
    ## abline(h=truePsi, col=1, lty=3)
    ## abline(v=truePsi, col=1, lty=3)
    abline(a=0, b=1, lty=2, col=1)
    ## usr <- par("usr")
    par(usr=c(usr[1:2], 0, 5), new=TRUE)
    plot(d0, axes=FALSE, xlab="", ylab="", main="", col=1, lty=1, xlim=xlim)
    par(usr=c(usr[1:2], 0, 5), new=TRUE)
    plot(d1$y, d1$x, axes=FALSE, xlab="", ylab="", main="", type="l", col=2, lty=1, ylim=xlim)
    if (doGraphs) {
      dev.off()
    }
    ##
    div <- sapply(NPVI, function(mat){mat[, "div"]})

    ##
    ## misc.
    ##
    (l1[steps][1]-l1[steps][-1])/l1[steps][1]
    
    ## Checking normality of estimate
    ##
    ## ## with theoretical putative asymptotic variance 
    ksPs <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-truePsi)/sqrt(sim0$varIC), pnorm)$p.value)
    ksPs.lower <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-confInt0[1])/sqrt(sim0$varIC), pnorm)$p.value) ## using lower end of confidence interval on truePsi 
    ksPs.upper <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-confInt0[2])/sqrt(sim0$varIC), pnorm)$p.value) ## using upper end of confidence interval on truePsi 
    ## ## with estimated putative asymptotic variance                                                                         
    ksPs2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-truePsi)/sic[ii, ], pnorm)$p.value)
    ksPs.lower2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-confInt0[1])/sic[ii, ], pnorm)$p.value) ## using lower end of confidence interval on truePsi 
    ksPs.upper2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-confInt0[2])/sic[ii, ], pnorm)$p.value) ## using upper end of confidence interval on truePsi 

    plot(ksPs, ylim=c(0, max(ksPs)))
    
    ksPs
    ksP <- ksPs[steps]
    ksP.lower <- ksPs.lower[steps]
    ksP.upper <- ksPs.upper[steps]
    ksP2 <- ksPs2[steps]
    ksP.lower2 <- ksPs.lower2[steps]
    ksP.upper2 <- ksPs.upper2[steps]
    ## all very small: retrospectively, due to remainding bias
    
    ## Retrieving the test statistics
    ## ## with theoretical putative asymptotic variance 
    ksStats <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-truePsi)/sqrt(sim0$varIC), pnorm)$statistic)
    ksStats.lower <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-confInt0[1])/sqrt(sim0$varIC), pnorm)$statistic) ## using lower end of confidence interval on truePsi 
    ksStats.upper <- apply(psi, 1, FUN=function(x) ks.test(sqrt(nObs)*(x-confInt0[2])/sqrt(sim0$varIC), pnorm)$statistic) ## using upper end of confidence interval on truePsi 
    ## ## with estimated putative asymptotic variance                                                                         
    ksStats2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-truePsi)/sic[ii, ], pnorm)$statistic)
    ksStats.lower2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-confInt0[1])/sic[ii, ], pnorm)$statistic) ## using lower end of confidence interval on truePsi 
    ksStats.upper2 <- sapply(1:nrow(psi), FUN=function(ii) ks.test(sqrt(nObs)*(psi[ii, ]-confInt0[2])/sic[ii, ], pnorm)$statistic) ## using upper end of confidence interval on truePsi 

    ksStat <- ksStats[steps]
    ksStat.lower <- ksStats.lower[steps]
    ksStat.upper <- ksStats.upper[steps]
    ksStat2 <- ksStats2[steps]
    ksStat.lower2 <- ksStats.lower2[steps]
    ksStat.upper2 <- ksStats.upper2[steps]

    ksZ <- qnorm(1-ksP/2)
    
    ## Checking coverage

    ## not good: this is not empirical coverage as it relies on sim0 !
    confInt <- truePsi+c(-1, 1)*qnorm(.975)*sqrt(sim0$varIC/nObs)

    ## the good way
    tol <- qnorm(.975)*sic/sqrt(nObs)

    covered <- (psi > truePsi-tol) & (psi < truePsi+tol)
    covered.lower <- (psi > confInt0[1]-tol) & (psi < confInt0[1]+tol)
    covered.upper <- (psi > confInt0[2]-tol) & (psi < confInt0[2]+tol)

    coverage <- apply(covered, 1, mean)
    coverage.lower <- apply(covered.lower, 1, mean)
    coverage.upper <- apply(covered.upper, 1, mean)

    coverage.optimistic <- apply(covered | covered.lower | covered.upper, 1, mean)
    
    coverage <- round(coverage[steps], 3)
    coverage.lower <- round(coverage.lower[steps], 3)
    coverage.upper <- round(coverage.upper[steps], 3)
    coverage.optimistic <- round(coverage.optimistic[steps], 3)
    
    tab <- rbind(ksStat, ksP, coverage)
    nms <- paste('$\\psi_n^', 1:ncol(tab)-1, "$", sep="") ## does not really work
    rownames(tab) <- c("KS test statistic", "p value of KS test", "Empirical coverage (95%)")

    library(xtable)
    tabPath <- "tex"
    tabPath <- Arguments$getWritablePath(tabPath)
    filename <- sprintf("table,%s.tex", figTag)
    pathname <- file.path(tabPath, filename)
    xt <- xtable(tab, display=c("s", rep("g", length(steps))))
    print(xt, file=pathname)

  }
}

