estimateG <- function(obs, weights=NULL,
                      flavor=c("learning", "superLearning", "h2oEnsembleLearning"), learnG,
                      light=TRUE, SuperLearner.=NULL, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);

  ## Argument 'weights':
  weights <- validateArgumentObsWeights(weights, nrow(obs));
  
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

  ## Argument 'SuperLearner.'
  if (flavor!="learning") {
    if (is.null(SuperLearner.) || mode(SuperLearner.)!="function") {
      throw("Argument 'SuperLearner.' should be a function")
    }
  }

  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  if (flavor=="learning") {
    g <- learnG(obs, weights=weights, light=light, ...);
  } else if (flavor=="superLearning") {
    obsD <- as.data.frame(obs)
    logSL <- as.logical(less(verbose, 10));  ## decrease verbosity in SuperLearner
    SL.library.g <- learnG;

    fitG <- SuperLearner.(Y=(obsD[, "X"]==0)+0, X=extractW(obsD), ## obsD[, "W", drop=FALSE]
                          obsWeights=weights,
                          SL.library=SL.library.g, verbose=logSL,
                          family=binomial(), ...)
    g <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitG, newdata=Wd)$pred
    }  
  } else if (flavor=="h2oEnsembleLearning") {
    EL.library.g <- learnG;
    obsD <- as.data.frame(obs)
    obsD$Y <- as.factor(as.integer(obsD[, "X"]==0)) ## forces binary classification
    data <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), obsD)
    
    fitG <- SuperLearner.(y="Y", x=colnames(extractW(obsD)),
                          training_frame=data,
                          family="binomial",
                          learner=EL.library.g,
                          weights_column=weights)
    g <- function(W) {
      Wd <- as.data.frame(W)
      newdata <- h2o::as.h2o(attr(SuperLearner., "H2OConnection"), Wd)
      predict(fitG, newdata=newdata)$pred
    }  
  }

  verbose && cat(verbose, "g(W):");
  verbose && print(verbose, summary(g(extractW(obs))));

  g
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

