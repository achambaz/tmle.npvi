setMethodS3("updatePsi", "NPVI", function(this, B, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'B':
  B <- Arguments$getInteger(B);
  
  ## Argument 'verbose'
  verbose <- Arguments$getVerbose(verbose);

  ## Retrieve parameters
  parsimonious <- getParsimonious(this)
  fX <- getFX(this)
  obs <- getObs(this)
  obsWeights <- getObsWeights(this)
  W <- obs[, "W"]
  X <- fX(obs)
  Xq <- getXq(this)
  
  g <- getG(this);
  mu <- getMu(this);
  muAux <- getMuAux(this);
  sigma2 <- getSigma2(this);

  obsWeights <- getObsWeights(this);
  weightsW <- getWeightsW(this);

  ## Perform 'B' simulations according to the estimated parameters
  verbose && enter(verbose, "Simulating ", B, " observations");
  obsB <- simulateData(B, W, X, Xq, g, mu, muAux, sigma2,
                       weights=obsWeights,
                       weightsW=weightsW, parsimonious=parsimonious, verbose=verbose)
  verbose && str(verbose, obsB);
  verbose && exit(verbose);

  ## Compute 'theta' and 'theta0' on these B samples
  theta <- getTheta(this)
  theta0 <- getTheta0(this)

  ## Estimate psiPn:
  psi1 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obs,
                      sigma2=sigma2, weights=obsWeights, verbose=verbose) 
  this$.psiPn <- psi1$mean;
  this$.psiPn.sd <- psi1$sd;
  ## Estimate psi:
  psi0 <- estimatePsi(theta=theta, theta0=theta0, fX=fX, obs=obsB,
                      sigma2=sigma2, weights=rep(1/B, B), verbose=verbose) 
  this$.psi <- psi0$mean;
  this$.psi.sd <- psi0$sd;
### CAUTION
### CAUTION: interpretation of 'this$.psiPn.sd' dubious, contrary to that
###          of 'this$.psi.sd'
### CAUTION

  rm(obsB)
})

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

