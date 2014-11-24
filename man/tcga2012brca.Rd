\name{tcga2012brca}
\alias{tcga2012brca}
\docType{data}
\title{Sample breast cancer data from TCGA}
\description{
  Expression, DNA copy number, and DNA methylation data of 125
  genes of chromosome 21 for 463 breast cancer samples from TCGA.
}
\usage{tcga2012brca}
\format{
A list of matrices.  Each matrix corresponds to a particular gene and
has the following columns:\describe{
    \item{\eqn{Y}}{gene expression level (outcome  variable)}
    \item{\eqn{X}}{DNA copy number (continuous exposure variable)}
    \item{\eqn{W1 \dots Wk}}{DNA methylation (k baseline covariates)}
  }
}
\details{
  These data were obtained using the scripts located in
  \code{testScripts/tcga2012brca}. See the \code{tmle-npvi.pdf} vignette
  for more details on the preparation of the data set.

  Gene names and genomic coordinates are stored in the names of the list.
}
\source{
  Raw data may be retrieved from \url{https://tcga-data.nci.nih.gov/docs/publications/brca_2012/}

  \describe{
\item{methylation}{\url{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.methylation.27k.450k.466.zip}}
\item{copy number data}{\url{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.GISTIC2.tar.gz}}
\item{expression}{\url{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.exp.466.med.txt}}
\item{annotation data for methylation probes}{\url{
http://supportres.illumina.com/documents/myillumina/b78d361a-def5-4adb-ab38-e8990625f053/humanmethylation450_15017482_v1-2.csv}
}}}

\references{Chambaz, A., Neuvial, P., & van der Laan, M. J. (2012).
Estimation of a non-parametric variable importance measure of a
continuous exposure. Electronic journal of statistics, 6, 1059--1099.

Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. Nature, 490(7418), 61-70.}

\author{Antoine Chambaz, Pierre Neuvial}
\examples{
  data(tcga2012brca)
  nms <- names(tcga2012brca)
  ii <- grep("TP53", nms)
  obs <- tcga2012brca[[ii]]
  pairs(obs, main=nms[ii])

  thr <- 0.02
  whichSmall <- which(abs(obs[, "X"]) <= thr)
  obs[whichSmall, "X"] <- 0

  ## the code below takes ~20s which is above the acceptable limit for R CMD check
  \dontrun{res <- tmle.npvi(obs)}
}
\keyword{datasets}
