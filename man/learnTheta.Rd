\name{learnTheta}
\alias{learnTheta}
\title{Estimation of Cond. Expect. of Y given (X,W)}
\description{Function  for the  estimation of  the conditional  expectation  of \eqn{Y}
given \eqn{(X, W)}, version based on 'glm'.}
\arguments{
  \item{obs}{The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
argument of the  \code{function} \code{tmle.npvi}, where only observations
with \eqn{X\neq 0} are kept.}
  \item{light}{A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
\code{TRUE} (its default value). This requires that the result of each fit
be reduced  in size (for  a faster execution). Currently  implemented only
for flavor \code{learning}.}
  \item{\dots}{Additional arguments possibly needed.}
}

\value{Returns the fitted object.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{learnG, learnMuAux, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta}

