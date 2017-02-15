\name{getPsi.NPVI}
\alias{getPsi.NPVI}
\alias{getPsi}
\title{Returns Current Estimator}
\description{Returns the current value of the estimator.}
\usage{\method{getPsi}{NPVI}(this, ...)}
\arguments{
  \item{this}{An object of class \code{TMLE.NPVI}.}
  \item{\dots}{Not used.}
}

\value{Retrieves  the current value  of the estimator \eqn{\Psi(P_n^k)}  of the
parameter  of interest. Its computation involves  simulation of a large number of 
iid copies of \eqn{(X,W)} under \eqn{P_n^k}.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{tmle.npvi, getHistory, getPsiSd}

