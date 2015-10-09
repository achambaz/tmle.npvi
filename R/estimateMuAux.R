estimateMuAux <- function(obs, weights, id,
                          flavor=c("learning", "superLearning", "h2oEnsembleLearning"),
                          learnMuAux,
                          light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'weights':
  weights <- Arguments$getNumerics(weights);  

  ## Argument 'id':
  id <- Arguments$getCharacters(id);  

  ## Argument 'flavor':
  flavor <- match.arg(flavor);
  learnMode <- switch(flavor,
                      learning="function",
                      superLearning="character",
                      h2oEnsembleLearning="character");

  ## Argument 'learnMuAux'
  mode <- mode(learnMuAux);
  if (mode != learnMode) {
    throw("Argument 'learnMuAux' should be of mode '", learnMode, "', not '", mode, "' for flavor: ", flavor);
  }

  ## Argument 'SuperLearner.'
  if (flavor!="learning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  idx <- which(obs[, "X"] != 0);
  
  if (flavor=="learning") {
    muAux <- learnMuAux(obs[idx, ], weights[idx], light=light, ...);
  } else if (flavor=="superLearning") {
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.muAux <- learnMuAux;
    obsD <- as.data.frame(obs)
    ## WW <- obsD[idx, "W", drop=FALSE]
    WW <- extractW(obsD[idx, ])

    fitMuAux <- SuperLearner.(Y=obsD[idx, "X"], X=WW,
                              obsWeights=weights[idx], id=id,
                              SL.library=SL.library.muAux, verbose=logSL,
                              family=gaussian(), ...);
    verbose && print(verbose, fitMuAux);
    muAux <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitMuAux, newdata=Wd)$pred;
    }
  } else if (flavor=="h2oEnsembleLearning") {
    EL.library.muAux <- learnMuAux;
    obsD <- as.data.frame(obs)
    obsD <- obsD[idx, ]
    data <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), obsD)

    ##
    ## CAUTION: provide 'id' as soon as this argument is supported
    ##
    
    fitMuAux <- SuperLearner.(y="X", x=colnames(extractW(obsD)),
                              training_frame=data,
                              family="gaussian",
                              learner=EL.library.muAux,
                              weights_column=weights)
    muAux <- function(W) {
      Wd <- as.data.frame(W)
      newdata <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), Wd)
      predict(fitMuAux, newdata=newdata)$pred;
    }
  }
  verbose && cat(verbose, "mu'(W):");
  verbose && print(verbose, summary(muAux(extractW(obs))));
  
  muAux
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

