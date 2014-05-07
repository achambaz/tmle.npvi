library(tmle.npvi)
library("R.utils")
log <- Arguments$getVerbose(-8, timestamp=TRUE)

path <- "geneData/tcga_brca_2012"
path <- Arguments$getReadablePath(path)
files <- list.files(path)

TMLE <- vector("list", length(files))
names(TMLE) <- unlist(strsplit(files, split=".xdr"))
for (ii in 1:1) {# length(files)) {
  ## loading the data
  pathname <- file.path(path, files[ii])
  obs <- loadObject(pathname)
  ## thresholding copy number data
  whichSmall <- which(abs(obs[, "X"]) <= 2e-2)
  obs[whichSmall, "X"] <- 0
  ##
  tmle <- try(tmle.npvi(obs, f=identity, flavor="learning"))
  if (inherits(tmle, "try-error")) {
    TMLE[[ii]] <- NA
  } else {
    TMLE[[ii]] <- tmle
  }
}

