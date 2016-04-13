estimateDevMu <- function(muW, obs, weights, id,
                          eic1, flavor=c("learning", "superLearning", "h2oEnsembleLearning"),
                          learnDevMu,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'mu':
  muW <- Arguments$getNumerics(muW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);

  ## Argument 'weights':
  obsWeights <- Arguments$getNumerics(weights);  

  ## Argument 'id':
  id <- Arguments$getCharacters(id);  
  
  ## Argument 'eic1'
  eic1 <- Arguments$getNumerics(eic1);
  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                         learning="function",
                         superLearning="character",
                         h2oEnsembleLearning="character");

  ## Argument 'learnDevMu'
  mode <- mode(learnDevMu);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevMu' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devMu <- learnDevMu(obs, obsWeights, eic1, muW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    SL.library.devMu <- learnDevMu
    obsD <- as.data.frame(obs)
    ZdevMu <- (obsD[, "X"] - muW) * eic1;

    fitDevMu <- SuperLearner.(Y=ZdevMu, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                              obsWeights=obsWeights, id=id,
                              SL.library=SL.library.devMu, verbose=logSL,
                              family=gaussian(), ...);
    devMu <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitDevMu, newdata=Wd)$pred;
    }
  } else if (flavor=="h2oEnsembleLearning") {
    EL.library.devMu <- learnDevMu;
    obsD <- as.data.frame(obs)
    ZdevMu <- (obsD[, "X"] - muW) * eic1;
    obsD$Y <- ZdevMu
    data <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), obsD)

    ##
    ## CAUTION: provide 'id' as soon as this argument is supported
    ##
    
    fitDevMu <- SuperLearner.(y="Y", x=colnames(extractW(obsD)),
                              training_frame=data,
                              family="gaussian",
                              learner=EL.library.devMu,
                              weights_column=obsWeights)

    devMu <- function(W) {
      Wd <- as.data.frame(W)
      newdata <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), Wd)
      predict(fitDevMu, newdata=newdata)$pred;
    }
  }

  verbose && cat(verbose, "devMu(W):");
  verbose && print(verbose, summary(devMu(extractW(obs))));
  
  devMu
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

