###########################################################################/**
# @RdocFunction estimateDevG
# @alias estimateDevG
#
# @title "Estimates the direction in which parameter 'g' should be updated"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{g}{A @function that estimates the conditional probability of
#     neutral copy number (X==0) given the DNA methylation level: P(X=0|W)}
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
#   \item{learnDevG}{If \code{flavor=="learning"}, a function for learning
#     the direction in which parameter \var{g} should be updated. If
#     \code{flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{SuperLearner::SuperLearner} for learning the 
#     direction in which parameter \var{g} should be updated.}
#   \item{\dots}{Further arguments to be passed to 'learnDevG' for the
#     "learning" flavor, and to 'SuperLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  A @function, the estimated \var{devG}.
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
estimateDevG <- function(gW, obs, eic1, flavor=c("learning", "superLearning"), learnDevG,
                         light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'gW':
  gW <- Arguments$getNumerics(gW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=FALSE);

  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                      learning="function",
                      superLearning="character");

  ## Argument 'learnDevG'
  mode <- mode(learnDevG);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevG' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devG <- learnDevG(obs, eic1, gW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    obsD <- as.data.frame(obs)
    ZdevG <- eic1 * ( (obsD[, "X"]==0) - gW );
    SL.library.devG <- learnDevG;

    fitDevG <- SuperLearner.(Y=ZdevG, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                             SL.library=SL.library.devG, verbose=logSL,
                             family=gaussian(), ...)
    devG <- function(W) {
      Wd <- as.data.frame(W)
      predict.SuperLearner(fitDevG, newdata=Wd)$pred
    }
  }
  verbose && cat(verbose, "devG(W):");
  verbose && print(verbose, summary(devG(extractW(obs))));

  devG
}


############################################################################
## HISTORY:
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-22
## o Created.
############################################################################

