estimatePsi <- function(theta, theta0, fX, obs, sigma2, weights, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'theta':
  mode <- mode(theta);
  if (mode != "function") {
    throw("Argument 'theta' should be of mode 'function', not '", mode);
  }

  ## Argument 'theta0':
  mode <- mode(theta0);
  if (mode != "function") {
    throw("Argument 'theta0' should be of mode 'function', not '", mode);
  }

  ## Argument 'fX':
  mode <- mode(fX);
  if (mode != "function") {
    throw("Argument 'fX' should be of mode 'function', not '", mode);
  }
  
  ## Argument 'obs':
  obs <- validateArgumentObs(obs);
  
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);

  ## Argument 'weights':
  obsWeights <- Arguments$getNumerics(weights);
  
  T <- theta(obs[, c("X", "W")]);
  verbose && cat(verbose, "theta(X, W):");
  verbose && str(verbose, T);
  T0 <- theta0(obs[, "W", drop=FALSE]);
  verbose && cat(verbose, "theta0(W):");
  verbose && str(verbose, T0);

  argument <- fX(obs) * (T - T0);
  mean.psi1 <- sum(argument*obsWeights);

  ## not valid if 'obsWeights' does not sum up to one...
  ## var.psi1 <- sum((argument^2)*obsWeights) - mean.psi1^2;
  ## instead,
  ## sum_b (Z_b - sum_beta Z_beta*w_beta )^2*w_b
  ## = sum_b Z_b^2*w_b - 2*(sum_b Z_b*w_b) + (sum_b w_b)*(sum_b Z_b*w_b),
  ## hence
  var.psi1 <- sum((argument^2)*obsWeights) - (2-sum(obsWeights))*mean.psi1^2;
### CAUTION
### CAUTION: dubious interpretation of 'var.psi1' when it is computed
###          based on 'obs' and not 'obsB'
### CAUTION
  mean.psi1 <- mean.psi1/sigma2;
  sd.psi1 <- sqrt(var.psi1)/sigma2;
    
  list(mean=mean.psi1, sd=sd.psi1);
}


############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

