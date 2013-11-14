###########################################################################/**
# @RdocFunction estimateMuAux
# @alias estimateMuAux
#
# @title "Estimates parameter 'muAux' from the observations"
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
#   \item{learnMuAux}{If \code{flavor=="learning"}, a function for learning
#     parameter \var{mu'=mu/(1-g)}, where \var{mu=E(X|W)} and \var{P(X=0|W)}. If
#     \code{flavor=="superLearning"}, a library of learning functions to be
#     passed to \code{SuperLearner::SuperLearner} for learning parameter
#     \var{mu'=mu/(1-g)}.}
#   \item{\dots}{Further arguments to be passed to 'learnMuAux' for the
#     "learning" flavor, and to 'SuperLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  Returns the estimated \var{muAux}.
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
estimateMuAux <- function(obs, flavor=c("learning", "superLearning"), learnMuAux,
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

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  idx <- which(obs[, "X"] != 0);
  
  if (flavor=="learning") {
    muAux <- learnMuAux(obs[idx, ], light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.muAux <- learnMuAux;
    obsD <- as.data.frame(obs)
    ## WW <- obsD[idx, "W", drop=FALSE]
    WW <- extractW(obsD[idx, ])

    fitMuAux <- SuperLearner.(Y=obsD[idx, "X"], X=WW,
                              SL.library=SL.library.muAux, verbose=logSL,
                              family=gaussian(), ...);
    verbose && print(verbose, fitMuAux);
    muAux <- function(W) {
      Wd <- as.data.frame(W)
      predict.SuperLearner(fitMuAux, newdata=Wd)$pred;
    }
  }
  verbose && cat(verbose, "mu'(W):");
  verbose && print(verbose, summary(muAux(extractW(obs))));
  
  muAux
}

############################################################################
## HISTORY:
## 2011-09-23
## o Now estimating only 'muAux' (thus independently of 'g').
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-17
## o Created.
############################################################################

