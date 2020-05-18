\name{learnCondExpX2givenW}
\alias{learnCondExpX2givenW}
\title{Estimation of Cond. Expect. of X^2 Given W}
\description{Function for  the estimation of  the conditional expectation  of \eqn{X^2}
given \eqn{W} when \code{flavor} is set to "learning".}
\arguments{
  \item{obs}{The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
argument of the \code{function} \code{tmle.npvi}.}
  \item{light}{A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
\code{TRUE} (its default value). This requires that the result of each fit
be reduced  in size (for  a faster execution). Currently  implemented only
for flavor \code{learning}.}
}

\value{Returns the fitted object.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{learnG, learnMuAux, learnTheta, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta}

