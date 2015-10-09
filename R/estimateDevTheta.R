estimateDevTheta <- function(thetaXW, obs, weights, id,
                             flavor=c("learning", "superLearning", "h2oEnsembleLearning"),
                             learnDevTheta,
                             light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'thetaXW':
  thetaXW <- Arguments$getNumerics(thetaXW);
    
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=TRUE);

  ## Argument 'weights':
  weights <- Arguments$getNumerics(weights);  

  ## Argument 'id':
  id <- Arguments$getCharacters(id);  

  
  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnDevMode <- switch(flavor,
                         learning="function",
                         superLearning="character",
                         h2oEnsembleLearning="character");

  ## Argument 'learnDevTheta'
  mode <- mode(learnDevTheta);
  if (mode != learnDevMode) {
    throw("Argument 'learnDevTheta' should be of mode '", learnDevMode, "', not '", mode, "' for flavor: ", flavor);
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
    devTheta <- learnDevTheta(obs, weights, thetaXW, light=light, verbose=verbose);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in superLearner
    obsD <- as.data.frame(obs)
    ZdevTheta <- (obsD[, "Y"]-thetaXW)^2;
    SL.library.devTheta <- learnDevTheta;

    fitDevTheta <- SuperLearner.(Y=ZdevTheta, X=extractXW(obsD),  ## obsD[, c("X", "W")]
                                 obsWeights=weights, id=id,
                                 SL.library=SL.library.devTheta, verbose=logSL,
                                 family=gaussian(), ...);
    
    devTheta <- function(XW) {
      XWd <- as.data.frame(XW)
      predict(fitDevTheta, newdata=XWd)$pred
    }
  } else if (flavor=="h2oEnsembleLearning") {
    EL.library.devTheta <- learnDevTheta;
    obsD <- as.data.frame(obs)
    ZdevTheta <- (obsD[, "Y"]-thetaXW)^2;
    obsD$Y <- ZdevTheta
    data <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), obsD)

    ##
    ## CAUTION: provide 'id' as soon as this argument is supported
    ##
    
    fitDevTheta <- SuperLearner.(y="Y", x=colnames(extractXW(obsD)),
                                 training_frame=data,
                                 family="gaussian",
                                 learner=EL.library.devTheta,
                                 weights_column=weights)
    
    devTheta <- function(XW) {
      XWd <- as.data.frame(XW)
      newdata <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), XWd)
      predict(fitDevTheta, newdata=newdata)$pred
    }
  }

  verbose && cat(verbose, "devTheta(XW):");
  verbose && print(verbose, summary(devTheta(extractXW(obs))));

  devTheta
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

