###########################################################################/**
# @RdocFunction estimateDevMu
# @alias estimateDevMu
#
# @title "Estimates the direction in which parameter 'mu' should be updated"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mu}{A @function that estimates the conditional expectation of
#       the DNA copy number given the DNA methylation level: E(X|W)}
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#     \item{eic1}{A @numeric @vector, the first component of the estimated
#       efficient influence curve of parameter \var{psi}.}
#   \item{flavor}{A @character, the type of estimation to be performed.
#     Two flavors are supported: "learning" and "superLearning".}
#   \item{learnDevMu}{If \code{flavor=="learning"}, a function for learning
#     the direction in which parameter \var{mu} should be updated. If
#     \code{flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{superLearner::superLearner} for learning the 
#     direction in which parameter \var{mu} should be updated.}
#   \item{\dots}{Further arguments to be passed to 'learnDevMu' for the
#     "learning" flavor, and to 'superLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  A @function, the estimated \var{devMu}.
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
estimateDevMu <- function(muW, obs, eic1, flavor=c("learning", "superLearning"), learnDevMu,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'mu':
  muW <- Arguments$getNumerics(muW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);
  
  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devMu <- learnDevMu(obs, eic1, muW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    SL.library.devMu <- learnDevMu
    obsD <- as.data.frame(obs)
    ZdevMu <- (obsD[, "X"] - muW) * eic1;

    fitDevMu <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                              SL.library=SL.library.devMu, verbose=logSL,
                              family=gaussian(), ...);
    devMu <- function(W) {
      Wd <- as.data.frame(W)
      predict.SuperLearner(fitDevMu, newdata=Wd)$pred;
    }
  }
  verbose && cat(verbose, "devMu(W):");
  verbose && print(verbose, summary(devMu(extractW(obs))));
  
  devMu
}


############################################################################
## HISTORY:
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-22
## o Created.
############################################################################

