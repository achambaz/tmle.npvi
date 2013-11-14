###########################################################################/**
# @set "class=NPVI"
# @RdocMethod update
#
# @title "Returns an updated estimate using Targeted Maximum Likelikood
#   Estimation (TMLE)"
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
#   \item{learnDevG{If \code{\flavor=="learning"}, a function for learning
#     the direction in which parameter \var{g} should be updated. If
#     \code{\flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{SuperLearner::SuperLearner} for learning the 
#     direction in which parameter \var{g} should be updated.}
#   \item{learnDevMu}{If \code{\flavor=="learning"}, a function for learning
#     the direction in which parameter \var{mu} should be updated. If
#     \code{\flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{SuperLearner::SuperLearner} for learning the 
#     direction in which parameter \var{mu} should be updated.}
#   \item{learnDevTheta}{If \code{\flavor=="learning"}, a function for learning
#     the direction in which parameter \var{theta} should be updated. If
#     \code{\flavor=="superLearning"}, a library of learning functions to
#     be passed to \code{SuperLearner::SuperLearner} for learning the 
#     direction in which parameter \var{theta} should be updated.}
#   \item{bound}{A @numeric, upper bound for the magnitude of the change
#     in the estimated parameter. This argument is passed on to \code{updateEpsilon}}
#   \item{B}{A @numeric @value, the number of observations to be simulated.}
#   \item{light}{If @TRUE, minimal information is stored for each fitted
#     model.  Currently only implemented for flavor 'learning'.}
#   \item{cleverCovTheta}{If @TRUE, theta is updated using the corresponding
#     "clever covariate".}
#   \item{\dots}{Arguments to be passed to 'estimateDev*' functions.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns a NPVI @object, containing the updated estimate of parameter \var{psi}.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @see "SuperLearner::SuperLearner"
#   @see "initializeEstimation"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("update", "NPVI", function(object,
                                       flavor=c("learning", "superLearning"),
                                       cvControl=NULL,
                                       learnDevG=NULL,
                                       learnDevMu=NULL,
                                       learnDevTheta=NULL,
                                       bound=1e-1, B=1e4,
                                       light=TRUE, 
                                       cleverCovTheta=TRUE,
                                       exact=TRUE, trueGMu=NULL, ...,
                                       verbose=FALSE) {

  this <- object; ## to please R CMD CHECK 
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'flavor'
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                         learning="function",
                         superLearning="character");
  
  ## Argument 'learnDevG'
  mode <- mode(learnDevG);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevG' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'learnDevTheta'
  mode <- mode(learnDevTheta);
  if ((mode != learnDevMode) & (!cleverCovTheta)) {
    throw("Argument 'learnDevTheta' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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

  ## Argument 'cleverCovTheta'
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact'
  exact <- Arguments$getLogical(exact);

  ## Argument 'trueGMu'
  useTrueGMu <- (!is.null(trueGMu))
  
  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  verbose <- less(verbose, 10);

  ## Incrementing the current value of 'step'
  this$.step <- getStep(this)+1;

  verbose && cat(verbose, "Iteration ", getStep(this), "\n");

  if (flavor=="superLearning") {
    test <- exists("SuperLearner.", mode="function", envir=.GlobalEnv)
    if (!test) {
      if (is.null(cvControl)) {
        cvControl <- SuperLearner.CV.control(V=10L)
      }
      SuperLearner. <- function(...) {
        warning("Setting 'V=10' in 'SuperLearner.'")
        SuperLearner(cvControl=SuperLearner.CV.control(V=10), ...)
      }
      assign("SuperLearner.", SuperLearner., envir=.GlobalEnv)
    }
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Retrieve elements of 'this'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obs <- getObs(this);
  obsT <- getObs(this, tabulate=FALSE)
  family <- getFamily(this);
  tabulate <- getTabulate(this)
  g <- getG(this);
  mu <- getMu(this);
  theta <- getTheta(this);
  theta0 <- getTheta0(this);
  sigma2 <- getSigma2(this);
  psi <- getPsi(this);
  fW <- getFW(this)
  fX <- getFX(this)
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update divergence (part 1/2)
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cleverCovTheta) {
    div <- NA  ## 'div' cannot be calculated without further (otherwise unnecessary) assumptions
  } else {
    weightsW <- getWeightsW(this)
    fY <- getFY(this)
    obsB <- simulateData(B, obs[, "W"], obsT[, "X"], g, mu, sigma2,
                         theta=theta, Y=obsT[, "Y"], weightsW=weightsW, family=family)
    ## taken from 'updateEfficientInfluenceCurve'
    thetaXW <- theta(obsB[, c("X", "W")]);
    theta0W <- theta0(obsB[, "W", drop=FALSE]);
    muW <- mu(obsB[, "W", drop=FALSE])
    gW <- g(obsB[, "W", drop = FALSE])
    
    X <- fX(obsB)
    Y <- fY(obsB)
          
    D1 <- X * (thetaXW - theta0W - X * psi);
    D2 <- (Y - thetaXW) * (X - muW/gW*(X==0));
    verbose && summary(verbose, D1);
    verbose && summary(verbose, D2);
    verbose && print(verbose, sigma2);
    
    eic1 <- D1 / sigma2;
    eic2 <- D2 / sigma2;
    eic <- eic1 + eic2;
    verbose && summary(verbose, eic);

    partialDiv <- mean(abs(eic))
    rm(eic1, eic2, eic);
  }
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Efficient influence curve
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  eic <- getEfficientInfluenceCurve(this);
  eic1 <- eic[, "eic1"];
  rm(eic);

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update the estimation of relevant components of the distribution
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the estimation of relevant components of the distribution");

  if (cleverCovTheta) {
    updateEpsilonTheta(this);
    ## Update 'theta', then update EIC and 'epsilon' accordinlgy
    updateTheta(this, NULL, cleverCovTheta=cleverCovTheta, exact=exact);
    updatePsi(this, B, verbose=verbose);
    updateEfficientInfluenceCurve(this, obs);
    
    theta <- getTheta(this)
    theta0 <- getTheta0(this)
    psi <- getPsi(this)
    eic <- getEfficientInfluenceCurve(this)
    eic1 <- eic[, "eic1"]
    rm(eic)

    ## Update history
    setDivergence(this, div)
    updateHistory(this);
    if (!useTrueGMu) {
      ## browser(skipCalls=2)
      updateEpsilon(this, cleverCovTheta=TRUE, bound=bound);  ## preparing update of 'mu' and 'g'
    }
  } else {  ## update 'theta' without using a clever covariate, and don't update EIC and 'epsilon'
    updateEpsilon(this, cleverCovTheta=FALSE, bound=bound);
    div <- abs(getEpsilon(this))*partialDiv
    thetaXW <- theta(obs[, c("X", "W")])
    devTheta <- estimateDevTheta(thetaXW, obsT, flavor=flavor, learnDevTheta=learnDevTheta, light=light, ..., verbose=verbose);
    updateTheta(this, devTheta, cleverCovTheta=cleverCovTheta, exact=exact);
  }

  effICW <- NULL;

  if (!useTrueGMu) {
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Update estimation of 'g' and 'mu'
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    ## To please R CMD CHECK
    learnCondExpX2givenW <- NULL; rm(learnCondExpX2givenW)
    learnCondExpXYgivenW <- NULL; rm(learnCondExpXYgivenW)

    if (flavor=="learning") {
      condExpX2givenW <- learnCondExpX2givenW(obsT, light=light); ## a 'true' function
      condExpXYgivenW <- learnCondExpXYgivenW(obsT, light=light); ## a 'true' function
    } else if (flavor=="superLearning") {
      logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
      SL.library.condExpX2givenW <- learnCondExpX2givenW; 
      SL.library.condExpXYgivenW <- learnCondExpXYgivenW; 
      obsD <- as.data.frame(obsT)
      fitCondExpX2givenW <- SuperLearner.(Y=obsD[, "X"]^2, X=extractW(obsD), ## obsD[, "W", drop=FALSE],
                                          SL.library=SL.library.condExpX2givenW, verbose=logSL,
                                          family=gaussian(), ...)
      fitCondExpXYgivenW <- SuperLearner.(Y=obsD[, "Y"]*obsD[, "X"], X=extractW(obsD), ## obsD[, "W", drop=FALSE],
                                          SL.library=SL.library.condExpXYgivenW, verbose=logSL,
                                          family=gaussian(), ...)
      condExpX2givenW <- function(W) {
        Wd <- as.data.frame(W)
        predict.SuperLearner(fitCondExpX2givenW, newdata=Wd)$pred
      }
      condExpXYgivenW <- function(W) {
        Wd <- as.data.frame(W)
        predict.SuperLearner(fitCondExpXYgivenW, newdata=Wd)$pred
      }
      verbose && cat(verbose, "E(X^2|W):");
      verbose && print(verbose, summary(condExpX2givenW(extractW(obsD))));
      ##
      verbose && cat(verbose, "E(XY|W):");
      verbose && print(verbose, summary(condExpXYgivenW(extractW(obsD))));
    }
    
    ## pasted from 'estimateEpsilon'
    theMin <- min(eic1)
    theMax <- max(eic1)
    if (theMin > 0) {
      interval <- c(-1/(1.001*theMax), 1e3)
    }
    if (theMax < 0) {
      interval <- c(-1e3, -1/(1.001*theMin))
    }
    if (theMin<=0 & theMax>=0) {
      interval <- c(-1/(1.001*theMax), -1/(1.001*theMin))
    }
   
    effICW <- function(W) {
      realW <- fW(cbind(W=W, X=NA, Y=NA))
      out <- (condExpXYgivenW(realW)-theta0(W)*mu(W)-psi*condExpX2givenW(realW))/sigma2;
      threshold(out, theMin, theMax)
    }
    
    
    ## Update 'mu' *before* 'g' as 'mu' depends on (the existing) 'g'.
    muW <- mu(extractW(obs))
    devMu <- estimateDevMu(muW, obsT, eic1, flavor=flavor, learnDevMu=learnDevMu, light=light, ..., verbose=verbose);
    updateMu(this, devMu, exact=exact, effICW=effICW);

    gW <- g(extractW(obs))
    devG <- estimateDevG(gW, obsT, eic1, flavor=flavor, learnDevG=learnDevG, light=light, ..., verbose=verbose);
    updateG(this, devG, exact=exact, effICW=effICW);

    ## Update 'sigma2'
    X <- fX(obs)
    devSigma2 <- mean(eic1 * X^2);
    updateSigma2(this, devSigma2);  
    verbose && exit(verbose);
    updateWeightsW(this, effICW)    
  }

  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update estimation of 'psi'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating the estimation of 'psi'");
  updatePsi(this, B, verbose=verbose);
  psi <- getPsi(this);
  verbose && cat(verbose, "psi: ", round(psi, 3));

  updateEfficientInfluenceCurve(this);  
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Update history
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setDivergence(this, div)
  updateHistory(this);

  verbose && exit(verbose);
})

############################################################################
## HISTORY:
## 2011-05-18
## o BUG FIX: the updated versions of 'theta', 'theta0' and 'psi' were not
## taken into account for the estimation of weights when cleverCov=TRUE.
## o 'effICW' is now thresholded to make sure that all weights are >0.
## 2011-02-28
## o Added update of 'step' element of structure.
## 2011-02-23
## o CLEANUP: rewrote NPVI.update using side functions 'estimateDev*' to
## estimate 'gradients'.
## 2011-02-15
## o Renamed EstimateECM into NPVI.
## 2011-02-08
## o Added parameters 'gmax', 'mumin', 'mumax', 'thetamin', 'thetamax'.
## 2011-02-07
## o Added parameter 'family' to list of arguments.
## 2011-02-04
## o Added 'psiPn'
## o Added option 'useTrueGMu'
## 2011-01-31
## o Added 'logLikIncr'
## 2011-01-24
## o Added argument 'exact' in order to exploit the "exact" version of the
##   updating process.
## 2011-01-20
## o Added argument 'cleverCovTheta' to use a "clever covariate" for updating
##   theta.
## 2010-12-31
## o Added argument 'gmin' for thresholding
## 2010-12-02
## o Added argument 'light' to store minimal information for fitted objects.
## 2010-11-29
## o Gradient learning now based only on the relevant part of the efficient
## influence curve.
## o Now using 'simulateData' to estimate psi.
## 2010-11-24
## o theta0 is now updated.
## 2010-05-01
## o Created.
############################################################################

