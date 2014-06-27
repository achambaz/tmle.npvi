## Parsimonious conditional simulation of X given W failed
nm <- "chr10,001095,WDR37"


library("tmle.npvi")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

path <- "geneData/tcga2012brca"
path <- Arguments$getReadablePath(path)
files <- list.files(path)

if (FALSE) {
  nas <- sapply(files, function(ff) {
    obs <- loadObject(file.path(path, ff))
    sum(is.na(obs))
  })
  sum(nas>0)
}

descr <- list(thresh=2e-2,
              f=identity,
              flavor="learning",
              iter=10,
              stoppingCriteria=list(mic = 0.001, div = 0.001, psi = 0.01))

idx <- grep(nm, files)
pathname <- file.path(path, files[idx])
obs <- loadObject(pathname)
nbcov <- ncol(extractW(obs))
if (nbcov==1) {
  colnames(obs) <- c("Y", "X", "W")
}

## thresholding copy number data
whichSmall <- which(abs(obs[, "X"]) <= descr$thresh)
obs[whichSmall, "X"] <- 0

##
tmle <- try(tmle.npvi(obs, f=descr$f, flavor=descr$flavor,
                      stoppingCriteria=descr$stoppingCriteria))

## mumax <- quantile(obs[, "X"], 0.9, type=1)
## mumin <- quantile(obs[, "X"], 0.1, type=1)
## tmle2 <- try(tmle.npvi(obs, f=descr$f, flavor=descr$flavor,
##                       stoppingCriteria=descr$stoppingCriteria, mumax=mumax, mumin=mumin))
