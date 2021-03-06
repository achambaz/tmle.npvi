setMethodS3("init", "NPVI", function(this, flavor=c("learning", "superLearning", "h2oEnsembleLearning"),
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
                      superLearning="character",
                      h2oEnsembleLearning="character");
  
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

  ## Retrieving 'obs'
  obs <- getObs(this, tabulate=FALSE);

  ## Retrieving 'weights'
  obsWeights <- getObsWeights(this);

  ## Retrieving 'id'
  id <- getId(this);
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## learning
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  verbose && enter(verbose, "Estimating relevant features of the distribution");
  
  if (!useTrueGMu) {
    g <- estimateG(obs, weights=obsWeights, id=id,
                   flavor=flavor, learnG=learnG, light=light,
                   SuperLearner.=SuperLearner.,
                   ..., verbose=verbose);
    muAux <- estimateMuAux(obs, weights=obsWeights, id=id,
                           flavor=flavor, learnMuAux=learnMuAux, light=light,
                           SuperLearner.=SuperLearner.,
                           ..., verbose=verbose);
  } else {
    g <- trueG
    muAux <- trueMuAux
  } 
  initializeG(this, g);
  initializeMu(this, muAux, g);

  theta <- estimateTheta(obs, weights=obsWeights, id=id,
                         flavor=flavor, learnTheta=learnTheta, light=light,
                         SuperLearner.=SuperLearner.,
                         ..., verbose=verbose);
  initializeTheta(this, theta);

  sigma2 <- sum(obs[, "X"]^2 * obsWeights);
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
## 2014-02-07
## o Created.
############################################################################

