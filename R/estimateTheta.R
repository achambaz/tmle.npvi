###########################################################################/**
# @RdocMethod estimateTheta
# @alias estimateTheta
#
# @title "Estimates parameter 'theta' from the observations"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{flavor}{A @character, the type of estimation to be performed.
#     Two flavors are supported: "learning" and "superLearning".}
#   \item{learnTheta}{If \code{\flavor=="learning"}, a function for learning
#     parameter \var{theta(X,W)=E[Y|X,W]}. If \code{\flavor=="superLearning"}, a library of
#     learning functions to be passed to \code{SuperLearner::SuperLearner}
#     for learning parameter \var{theta}.}
#   \item{\dots}{Further arguments to be passed to 'learnTheta' for the
#     "learning" flavor, and to 'SuperLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  Returns the estimated \var{theta}.
# }
#
# @author
#
# \seealso{
#   @seemethod "initialize.NPVI"
#   @seeclass "NPVI"
# }
#
#*/###########################################################################
estimateTheta <- function(obs, flavor=c("learning", "superLearning"), learnTheta,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnTheta'
  mode <- mode(learnTheta);
  if (mode != learnMode) {
    throw("Argument 'learnTheta' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  if (flavor=="learning") {
    theta <- learnTheta(obs, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.theta <- learnTheta;
    obsD <- as.data.frame(obs)

    fitTheta <- SuperLearner.(Y=obsD[, "Y"], X=extractXW(obsD), ## obsD[, c("X", "W")]
                              SL.library=SL.library.theta, verbose=logSL,
                              family=gaussian(), ...)
    theta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict.SuperLearner(fitTheta, newdata=XWd)$pred
    }
  }
  verbose && cat(verbose, "theta(X,W):");
  verbose && print(verbose, summary(theta(extractXW(obs)))); ## obs[, c("X", "W")]

  theta
}


############################################################################
## HISTORY:
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-17
## o Created.
############################################################################

