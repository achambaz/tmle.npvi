estimateDevG <- function(gW, obs, weights, id,
                         eic1, flavor=c("learning", "superLearning", "h2oEnsembleLearning"),
                         learnDevG,
                         light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'gW':
  gW <- Arguments$getNumerics(gW);
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=FALSE);

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

  ## Argument 'learnDevG'
  mode <- mode(learnDevG);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevG' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devG <- learnDevG(obs, obsWeights, eic1, gW, light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    obsD <- as.data.frame(obs)
    ZdevG <- eic1 * ( (obsD[, "X"]==0) - gW );
    SL.library.devG <- learnDevG;

    fitDevG <- SuperLearner.(Y=ZdevG, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                             obsWeights=obsWeights, id=id,
                             SL.library=SL.library.devG, verbose=logSL,
                             family=gaussian(), ...)
    devG <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitDevG, newdata=Wd)$pred
    }
  } else if (flavor=="h2oEnsembleLearning") {
    EL.library.devG <- learnDevG;
    obsD <- as.data.frame(obs)
    ZdevG <- eic1 * ( (obsD[, "X"]==0) - gW );
    obsD$Y <- ZdevG
    data <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), obsD)

    ##
    ## CAUTION: provide 'id' as soon as this argument is supported
    ##
    
    fitDevG <- SuperLearner.(y="Y", x=colnames(extractW(obsD)),
                             training_frame=data,
                             family="gaussian",
                             learner=EL.library.devG,
                             weights_column=obsWeights)

    devG <- function(W) {
      Wd <- as.data.frame(W)
      newdata <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), Wd)
      predict(fitDevG, newdata=newdata)$pred
    }
  }
  verbose && cat(verbose, "devG(W):");
  verbose && print(verbose, summary(devG(extractW(obs))));

  devG
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

