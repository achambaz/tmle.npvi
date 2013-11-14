###########################################################################/**
# @set "class=NPVI"
# @RdocMethod init
#
# @title "Performs an initial estimation of the parameter."
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{flavor}{A @character, the type of estimation to be performed.
#     Two flavors are supported: "learning" and "superLearning".}
#   \item{learnG}{If \code{flavor=="learning"}, a function for learning
#     parameter \var{g=P(X=0|W)}. If \code{flavor=="superLearning"}, a
#     library of learning functions to be passed to
#     \code{SuperLearner::SuperLearner} for learning parameter \var{g}.}
#   \item{learnMuAux}{If \code{\flavor=="learning"}, a function for learning
#     parameter \var{mu'=mu/(1-g)}, where \var{mu=E(X|W)}. If
#     \code{\flavor=="superLearning"}, a library of learning functions to be
#     passed to \code{SuperLearner::SuperLearner} for learning parameter
#     \var{mu'=mu/(1-g)}.}
#   \item{learnTheta}{If \code{\flavor=="learning"}, a function for learning
#     parameter \var{theta(X,W)=E[Y|X,W]}. If \code{\flavor=="superLearning"}, a library of
#     learning functions to be passed to \code{SuperLearner::SuperLearner}
#     for learning parameter \var{theta}.}
#   \item{bound}{A @numeric, upper bound for the magnitude of the change
#     in the estimated parameter. This argument is passed on to \code{updateEpsilon}}
#   \item{B}{A @numeric @value, the number of observations to be simulated.}
#   \item{family}{A @character, describing how X should be conditionally simulated
#     given W. Currently two options: 'parsimonious' and 'gaussian'.}
#   \item{light}{If @TRUE, minimal information is stored for each fitted
#     model.  Currently only implemented for flavor 'learning'.}
#   \item{\dots}{Arguments to be passed to 'estimate*' functions.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns an NPVI object containing estimates of \var{g}, \var{mu}, \var{theta},
#  and a corresponding estimate of parameter \var{psi}.
# }
#
# @author
#
# \seealso{
#   @see "NPVI.update"
#   @see "SuperLearner::SuperLearner"
#   @see "NPVI"
# }
#
#*/###########################################################################
setMethodS3("init", "NPVI", function(this, flavor=c("learning", "superLearning"),
                                     cvControl=NULL,
                                     learnG=NULL,
                                     learnMuAux=NULL,
                                     learnTheta=NULL,
                                     bound=1e-1, B=1e4,
                                     light=TRUE, 
                                     trueGMu=NULL,
                                     SuperLearner.=NULL,
                                     ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnTheta'
  mode <- mode(learnTheta);
  if (mode != learnMode) {
    throw("Argument 'learnTheta' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'bound':
  bound <- Arguments$getNumeric(bound);
  if (bound<=0) {
    throw("Argument 'bound' must be positive!\n")
  }

  ## Argument 'B':
  B <- Arguments$getInteger(B);
  
  ## Argument 'light'
  light <- Arguments$getLogical(light);

  ## Argument 'trueGMu'
  useTrueGMu <- (!is.null(trueGMu))
  if (useTrueGMu) {
    if (!is.list(trueGMu)) {
      throw("If not NULL, Argument 'trueGMu' should be a list")
    }
    trueG <- trueGMu[["g"]]
    if (mode(trueG) != "function") {
      throw("Argument 'trueGMu$g' should be a function, not a ", mode(trueG))
    }
    trueMuAux <- trueGMu[["muAux"]]
    if (mode(trueMuAux) != "function") {
      throw("Argument 'trueGMu$muAux' should be a function, not a ", mode(trueMuAux))
    }
  }

  ## Argument 'SuperLearner.'
  if (flavor=="superLearning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);
  verbose <- less(verbose, 10);

  ## retrieving 'obs'
  obs <- getObs(this, tabulate=FALSE);

  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## learning
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  verbose && enter(verbose, "Estimating relevant features of the distribution");
  
  if (!useTrueGMu) {
    g <- estimateG(obs, flavor=flavor, learnG=learnG, light=light,
                   SuperLearner.=SuperLearner.,
                   ..., verbose=verbose);
    muAux <- estimateMuAux(obs, flavor=flavor, learnMuAux=learnMuAux, light=light,
                           SuperLearner.=SuperLearner.,
                           ..., verbose=verbose);
  } else {
    g <- trueG
    muAux <- trueMuAux
  } 
  initializeG(this, g);
  initializeMu(this, muAux, g);

  theta <- estimateTheta(obs, flavor=flavor, learnTheta=learnTheta, light=light,
                         SuperLearner.=SuperLearner.,
                         ..., verbose=verbose);
  initializeTheta(this, theta);

  sigma2 <- mean(obs[, "X"]^2);
  setSigma2(this, sigma2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Updating 'psi' accordingly");
  updatePsi(this, B, verbose=verbose);
  psi0 <- getPsi(this);
  verbose && str(verbose, psi0);

  verbose && enter(verbose, "Updating efficient influence curve and 'epsilon' accordinlgy");
  updateEfficientInfluenceCurve(this);

  ## Update history
  updateHistory(this);

  verbose && exit(verbose);
})

############################################################################
## HISTORY:
## o Removed parameter 'family' and 'obs' from call to 'updatePsi'
##   and removed 'obs' and 'family' from arguments of 'init' method 
## o Converted function 'initializeEstimation' into a method
## 'initialize' for class NPVI.
## 2011-02-15
## o Renamed EstimateECM into NPVI.
## 2011-02-08
## o Added parameters 'gmax', 'mumin', 'mumax', 'thetamin', 'thetamax',
## the related functions, and adapted code accordingly.
## o Added thresholding of 'g', 'mu' and 'theta' IN "learning' MODE ONLY.
## 2011-02-07
## o Added parameter 'family' to arguments of function 'initializeEstimation'
##   and to attributes of its output 'res'
## 2011-02-04
## o Added parameter 'psiPn'
## o Added parameter 'useTrueGMu'
## 2011-01-24
## o Removed parameter 'theta0'.
## 2010-12-02
## o Added argument 'light' to store minimal information for fitted objects.
## 2010-11-29
## o Now using 'simulateData' to estimate psi.
## 2010-08-03
## o Added parameter 'theta0'.
## 2010-08-02
## o Added support for stantard (ie non-"super") learning.
## 2010-07-28
## o Fixed (known) bug in prediction of fitMuAux.
## 2010-07-09
## o Added a test to make sure that input data have a positive mass at 0.
## 2010-05-01
## o Created.
############################################################################

