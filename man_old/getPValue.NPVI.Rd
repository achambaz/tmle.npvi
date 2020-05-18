\name{getPValue.NPVI}
\alias{getPValue.NPVI}
\title{Calculates a p-value from a NPVI object}
\description{Calculates a p-value from a NPVI object}
\usage{\method{getPValue}{NPVI}(this, wrt.phi = TRUE, ...)}
\arguments{
  \item{this}{An object of class \code{TMLE.NPVI}.}
  \item{wrt.phi}{A  \code{logical}  equal  to  \code{TRUE}  by default,  which  means  that
\eqn{psi_n}  is  compared  with  \eqn{phi_n}.  Otherwise,  \eqn{psi_n}  is
compared with 0.}
  \item{\dots}{Not used.}
}

\value{Returns the p-value of the two-sided test of
``\eqn{Psi(P_0)=Phi(P_0)}'' or ``\eqn{Psi(P_0)=0}'', according to
the value of \code{wrt.phi}.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{tmle.npvi, getHistory, as.character.NPVI, getPValue.matrix}

