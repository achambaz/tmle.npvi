estimateTheta <- function(obs, weights, id,
                          flavor=c("learning", "superLearning"),
                          learnTheta,
                          familyY=c("gaussian", "binomial"),
                          light=TRUE, cvControl=NUlL, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);

  ## Argument 'weights':
  obsWeights <- Arguments$getNumerics(weights);  

  ## Argument 'id':
  id <- Arguments$getCharacters(id);  
  
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

  ## Argument 'familyY':
  familyY <- match.arg(familyY);

  ## Argument 'cvControl'
  if (flavor!="learning") {
    if (is.null(cvControl)) {
      throw("Argument 'cvControl' should have the same form as the output of 'SuperLearner.CV.control'")
    }
  }
  
  ## Argument 'SuperLearner.'
  if (flavor!="learning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  if (flavor=="learning") {
    theta <- learnTheta(obs, weights=obsWeights, light=light, family=familyY, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.theta <- learnTheta;
    obsD <- as.data.frame(obs)

    fitTheta <- SuperLearner.(Y=obsD[, "Y"], X=extractXW(obsD), ## obsD[, c("X", "W")]
                              obsWeights=obsWeights, id=id,
                              SL.library=SL.library.theta, verbose=logSL,
                              cvControl=cvControl,
                              family=familyY, ...)
    theta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict(fitTheta, newdata=XWd)$pred
    }
  }
  verbose && cat(verbose, "theta(X,W):");
  verbose && print(verbose, summary(theta(extractXW(obs)))); ## obs[, c("X", "W")]

  theta
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

