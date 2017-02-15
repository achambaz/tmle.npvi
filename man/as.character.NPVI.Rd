\name{as.character.NPVI}
\alias{as.character.NPVI}
\alias{as.character}
\title{Returns a Description}
\description{Returns a short string describing the NPVI object.}
\usage{\method{as.character}{NPVI}(x, ...)}
\arguments{
  \item{x}{An object of class \code{TMLE.NPVI}.}
  \item{\dots}{Not used.}
}

\value{A character string summarizing the content of the object. The summary contains:
\itemize{
\item{The sample size of the data set involved in the TMLE procedure.}
\item{The value of the TMLE and its estimated standard error.}
\item{A reminder of  the tuning of the stopping criteria,  and a report on
the convergence of the TMLE procedure (see \code{\link{tmle.npvi}}).}
\item{A confidence interval with default level of 95\% (the level can be changed by using \code{\link{setConfLevel}}).}
\item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=0}''.}
\item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=\Phi(P_0)}'', with the estimated value of \eqn{\Phi(P_0)}.}
}}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{tmle.npvi}

