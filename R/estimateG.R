###########################################################################/**
# @RdocFunction estimateG
# @alias estimateG
#
# @title "Estimates parameter 'g' from the observations"
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
#   \item{learnG}{If \code{flavor=="learning"}, a function for learning
#     parameter \var{g=P(X=0|W)}. If \code{flavor=="superLearning"}, a
#     library of learning functions to be passed to
#     \code{SuperLearner::SuperLearner} for learning parameter \var{g}.}
#   \item{\dots}{Further arguments to be passed to 'learnG' for the
#     "learning" flavor, and to 'SuperLearner' for the "superLearning"
#     flavor.}
# }
#
# \value{
#  The estimated \var{g}.
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
estimateG <- function(obs, flavor=c("learning", "superLearning"), learnG,
                      light=TRUE, ..., verbose=FALSE) {
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

  ## Argument 'learnG'
  mode <- mode(learnG);
  if (mode != learnMode) {
    throw("Argument 'learnG' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  if (flavor=="learning") {
    g <- learnG(obs, light=light, ...);
  } else if (flavor=="superLearning") {
    obsD <- as.data.frame(obs)
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.g <- learnG;

    ## To please R CMD CHECK
    SuperLearner. <- NULL; rm(SuperLearner.)
    fitG <- SuperLearner.(Y=(obsD[, "X"]==0)+0, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                          SL.library=SL.library.g, verbose=logSL,
                          family=binomial(), ...)
    g <- function(W) {
      Wd <- as.data.frame(W)
      predict.SuperLearner(fitG, newdata=Wd)$pred
    }  
  }
  verbose && cat(verbose, "g(W):");
  verbose && print(verbose, summary(g(extractW(obs))));

  g
}


############################################################################
## HISTORY:
## 2011-04-22
## o Call to new function 'validateArgumentObs'
## 2011-02-17
## o Created.
############################################################################

