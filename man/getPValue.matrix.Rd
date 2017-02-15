\name{getPValue.matrix}
\alias{getPValue.matrix}
\alias{getPValue}
\title{Calculates a p-value from a matrix object of type 'history'}
\description{Calculates a p-value from a matrix object of type 'history'}
\usage{\method{getPValue}{matrix}(this, wrt.phi = TRUE, nobs, ...)}
\arguments{
  \item{this}{The \code{history} of a TMLE procedure.}
  \item{wrt.phi}{A  \code{logical}  equal  to  \code{TRUE}  by default,  which  means  that
\eqn{psi_n}  is  compared  with  \eqn{phi_n}.  Otherwise,  \eqn{psi_n}  is
compared with 0.}
  \item{nobs}{An \code{integer}, the associated number of observations/clusters.}
  \item{\dots}{Not used.}
}

\value{Returns the p-value of the two-sided test of
``\eqn{Psi(P_0)=Phi(P_0)}'' of ``\eqn{Psi(P_0)=0}'', according to
the value of \code{wrt.phi}.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{tmle.npvi, getHistory.NPVI, as.character.NPVI}

