###########################################################################/**
# @RdocFunction estimateDevTheta
# @alias estimateDevTheta
#
# @title "Estimates the direction in which parameter 'theta' should be updated"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{theta}{A @function that estimates the conditional expectation of
#     the expression level given DNA copy number and DNA methylation: E(Y|X,W)}
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{flavor}{A @character, the type of estimation to be performed.
#     Two flavors are supported: "learning" and "superLearning".}
#   \item{learnDevTheta}{If \code{flavor=="learning"}, a function for learning
#     the direction in which parameter \var{theta} should be updated. If
#     \code{flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{superLearner::superLearner} for learning the 
#     direction in which parameter \var{theta} should be updated.}
#   \item{\dots}{Further arguments to be passed to 'learnDevTheta' for the
#     "learning" flavor, and to 'superLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  A @function, the estimated \var{devTheta}.
# }
#
# @author
#
# \seealso{
#   @seemethod "update.NPVI"
#   @seeclass "NPVI"
# }
#
#*/###########################################################################
estimateDevTheta <- function(thetaXW, obs, flavor=c("learning", "superLearning"), learnDevTheta,
                             light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'thetaXW':
  thetaXW <- Arguments$getNumerics(thetaXW);
    
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevTheta'
  mode <- mode(learnDevTheta);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevTheta' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devTheta <- learnDevTheta(obs, thetaXW, light=light, verbose=verbose);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    obsD <- as.data.frame(obs)
    ZdevTheta <- (obsD[, "Y"]-thetaXW)^2;
    SL.library.devTheta <- learnDevTheta;

    fitDevTheta <- SuperLearner.(Y=ZdevTheta, X=extractXW(obsD),  ## obsD[, c("X", "W")]
                                 SL.library=SL.library.devTheta, verbose=logSL,
                                 family=gaussian(), ...);
    
    devTheta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict.SuperLearner(fitDevTheta, newdata=XWd)$pred
    }
  }
  verbose && cat(verbose, "devTheta(XW):");
  verbose && print(verbose, summary(devTheta(extractXW(obs))));

  devTheta
}


############################################################################
## HISTORY:
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-22
## o Created.
############################################################################

