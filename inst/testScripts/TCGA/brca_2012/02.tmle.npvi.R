library(tmle.npvi)
library("R.utils")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

path <- "geneData/tcga_brca_2012"
path <- Arguments$getReadablePath(path)
files <- list.files(path)

nas <- sapply(files[idxs], function(ff) {
  obs <- loadObject(file.path(path, ff))
  sum(is.na(obs))
})

where <- unlist(strsplit(Sys.info()["nodename"], split="\\."))[1]
if (where!="ondine") {
  cArgs <- commandArgs()
  chunk <- as.character(cArgs[5])
  idx <- eval(parse(text=sub("-", ":", chunk)))
} else {
  idx <- 1:length(files)
  chunk <- paste(as.character(idx[c(1, length(idx))]), collapse="-")
}

files.idx <- files[idx]

descr <- list(thresh=2e-2,
              f=identity,
              flavor="learning",
              iter=10,
              stoppingCriteria=list(mic = 0.001, div = 0.001, psi = 0.01))

fileout <- paste(descr$flavor, "chromosome21", chunk, "RData", sep=".")

TMLE <- vector("list", length(files.idx))
names(TMLE) <- unlist(strsplit(files.idx, split=".xdr"))

counter <- 0
for (ii in 1:length(files.idx)) {
  counter <- counter+1
  ## loading the data
  pathname <- file.path(path, files.idx[ii])
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
  if (inherits(tmle, "try-error")) {
    TMLE[[ii]] <- attr(tmle, "condition")
  } else {
    TMLE[[ii]] <- list(nbcov=nbcov,
                       hist=getHistory(tmle))
  }
  ## saving
  if (counter==10) {
    save(descr, TMLE, file=fileout)
    counter <- 0
  }
}

save(descr, TMLE, file=fileout)
