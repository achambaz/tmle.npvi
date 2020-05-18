estimateG <- function(obs, weights, id,
                      flavor=c("learning", "superLearning"), learnG,
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
                      superLearning="character");

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
                          obsWeights=weights, id=id,
                          SL.library=SL.library.g, verbose=logSL,
                          family=binomial(), ...)
    g <- function(W) {
      Wd <- as.data.frame(W)
      predict(fitG, newdata=Wd)$pred
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

