validateArgumentObsWeights <- function(weights, nobs) {
  ## Argument 'nobs'
  nobs <- Arguments$getInteger(nobs, range=c(1, Inf))
  
  if (is.null(weights)) {
    weights <- rep(1/nobs, nobs);
  } else {
    weights <- Arguments$getNumerics(weights, range=c(0, 1));
    if (length(weights)!=nobs) {
      throw("Argument 'weights' should be a vector of size ", nobs, ", not ", length(weights));
    }
    sumWeights <- zapsmall(sum(weights))
    if (sumWeights!=1) {
      throw("Argument 'weights' should be a vector of non-negative numbers summing up to 1, not ", sumWeights); 
    }
  }
  weights;
}

############################################################################
## HISTORY:
## 2015-08-03
## o Created.
############################################################################

