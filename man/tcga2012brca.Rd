\name{tcga2012brca}
\alias{tcga2012brca}
\docType{data}
\title{Sample breast cancer data from TCGA}
\description{
  Expression, DNA copy number, and DNA methylation data of 10 genes for
  463 breast cancer samples from TCGA.
}
\usage{tcga2012brca}
\format{
A list of matrices.  Each matrix corresponds to a particular gene and
has the following columns:\description{
    \item{\eqn{Y}}{gene expression level (outcome  variable)}
    \item{\eqn{X}}{DNA copy number (continuous exposure variable)}
    \item{\eqn{W1 \dots Wk}}{DNA methylation (k baseline covariates)}
  }
}
\details{
  These data were obtained using the scripts located in
  \code{testScripts/tcga2012brca}. 
}
\source{
  Raw data may be retrieved from \url{https://tcga-data.nci.nih.gov/docs/publications/brca_2012/}

  \describe{
\item{methylation}{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.methylation.27k.450k.466.zip}
\item{copy number data}{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.GISTIC2.tar.gz}
\item{expression}{http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.exp.466.med.txt}
\item{annotation data for methylation probes}{
http://supportres.illumina.com/documents/myillumina/b78d361a-def5-4adb-ab38-e8990625f053/humanmethylation450_15017482_v1-2.csv
}}}

\references{Chambaz, A., Neuvial, P., & van der Laan, M. J. (2012).
Estimation of a non-parametric variable importance measure of a
continuous exposure. Electronic journal of statistics, 6, 1059--1099.

Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. Nature, 490(7418), 61-70.}

\author{Antoine Chambaz, Pierre Neuvial}
\examples{
  data(tcga2012brca)
  nms <- names(tcga2012brca)
  ii <- 1
  obs <- tcga2012brca[[ii]] 
  pairs(obs, main=nms[ii])
  res <- tmle.npvi(obs)
}
\keyword{datasets}
