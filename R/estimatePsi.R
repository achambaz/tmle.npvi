###########################################################################/**
# @RdocFunction estimatePsi
#
# @title "Returns an estimate of 'psi'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{theta}{A @function that estimates the conditional expectation of
#     the expression level given DNA copy number (centered) and DNA
#     methylation: E(Y|X,W).}
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{sigma2}{A @numeric, the squared expectation of DNA copy number
#     relative to the neutral state. If not provided, then
#     mean(X^2) is used.}
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns a @numeric, the estimated value \var{psi}.
# }
#
# @author
#
# \seealso{
#   @see "update"
#   @see "SuperLearner::SuperLearner"
#   @see "NPVI"
# }
#
#*/###########################################################################
estimatePsi <- function(theta, theta0, fX, obs, sigma2, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'theta':
  mode <- mode(theta);
  if (mode != "function") {
    throw("Argument 'theta' should be of mode 'function', not '", mode);
  }

  ## Argument 'theta0':
  mode <- mode(theta0);
  if (mode != "function") {
    throw("Argument 'theta0' should be of mode 'function', not '", mode);
  }

  ## Argument 'fX':
  mode <- mode(fX);
  if (mode != "function") {
    throw("Argument 'fX' should be of mode 'function', not '", mode);
  }
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);


  T <- theta(obs[, c("X", "W")]);
  verbose && cat(verbose, "theta(X, W):");
  verbose && str(verbose, T);
  T0 <- theta0(obs[, "W", drop=FALSE]);
  verbose && cat(verbose, "theta0(W):");
  verbose && str(verbose, T0);

  mean.psi1 <- mean(fX(obs) * (T - T0))/sigma2;
  sd.psi1 <- sd(fX(obs) * (T - T0))/sigma2;
    
  list(mean=mean.psi1, sd=sd.psi1);
}


############################################################################
## HISTORY:
## 2011-04-22
## o Removed parameter 'obsX' (sic!).
## o Added parameter 'fX'
## 2011-04-08
## o Added 'obsX' parameter.
## 2011-04-08
## o Added parameter 'theta0' (sic!).
## 2011-01-24
## o Removed parameter 'theta0'.
## 2010-11-26
## o Removed parameter 'X0'.  X is now assumed to be centered.
## 2010-08-03
## o Added parameter 'theta0'.
## 2010-07-12
## o If argument 'sigma2' is not provided, then mean(X^2) is used as an
##   estimate  of 'sigma2'.
## 2010-07-09
## o Added argument 'X0' to allow estimation of psi when copy number data
##   are not centered.
## 2010-05-01
## o Created.
############################################################################

