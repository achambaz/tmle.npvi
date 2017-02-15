\name{learnG}
\alias{learnG}
\title{Estimation of Cond. Prob. of X=x_0 Given W}
\description{Function  for   the  estimation   of  the  conditional   probability  that
\eqn{X=x_0} (the reference value for \eqn{X}) given \eqn{W}, version based
on 'glm'.}
\arguments{
  \item{obs}{The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
argument of the \code{function} \code{tmle.npvi}.}
  \item{theX0}{The reference value for \eqn{X}.}
  \item{light}{A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
\code{TRUE} (its default value). This requires that the result of each fit
be reduced  in size (for  a faster execution). Currently  implemented only
for flavor \code{learning}.}
  \item{\dots}{Additional arguments possibly needed. }
}

\value{Returns the fitted object.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta}

