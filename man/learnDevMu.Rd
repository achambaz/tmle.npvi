\name{learnDevMu}
\alias{learnDevMu}
\title{Estimation of Cond. Expect. of (X-muW)*effIC1 Given W }
\description{Function   for  the   estimation   of  the   conditional  expectation   of
\code{(X-muW)*effIC1}   given  \eqn{W}  when   \code{flavor}  is   set  to
"learning".}
\arguments{
  \item{obs}{The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
argument of the \code{function} \code{tmle.npvi}.}
  \item{effIC1}{The \code{vector}  \code{effIC1} of the  first component of  the efficient
influence curve, as currently estimated, evaluated at our observations.}
  \item{muW}{The  \code{vector} \code{muW}  of the  conditional expectation  of \eqn{X}
given \eqn{W}, as currently estimated, evaluated at our observations.}
  \item{light}{A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
\code{TRUE} (its default value). This requires that the result of each fit
be reduced  in size (for  a faster execution). Currently  implemented only
for flavor \code{learning}.}
  \item{verbose}{Prescribes the amount of information  output by the function.  Defaults to
\code{FALSE}.}
  \item{\dots}{Additional arguments possibly needed.}
}

\value{Returns the fitted object.}

\author{Antoine Chambaz, Pierre Neuvial}



\seealso{learnG, learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevTheta}

