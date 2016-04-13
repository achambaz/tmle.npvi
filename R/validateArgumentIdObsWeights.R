validateArgumentIdObsWeights <- function(id, weights, nobs) {
  ## Argument 'id'
  id <- Arguments$getCharacters(id);
  ## Argument 'nobs'
  nobs <- Arguments$getInteger(nobs, range=c(1, Inf));
  
  if (is.null(id)) {
    id <- as.character(1:nobs);
  } else {
    if (length(id)!=nobs) {
      throw("Argument 'id' should be a vector of size ", nobs, ", not ", length(id));
    }
  }
  
  if (is.null(weights)) {
    tab <- table(id);
    weights <- 1/tab[id];
    weights <- weights/sum(weights);
    attr(weights, "sum to one") <- TRUE
  } else {
    weights <- Arguments$getNumerics(weights, range=c(0, Inf));
    if (length(weights)!=nobs) {
      throw("Argument 'weights' should be a vector of size ", nobs, ", not ", length(weights));
    }
    sumWeights <- zapsmall(sum(weights));
    weights <- as.vector(weights);
    if (sumWeights!=1) {
      ## throw("Argument 'weights' should be a vector of non-negative numbers summing up to 1, not ", sumWeights);
      attr(weights, "sum to one") <- FALSE
    } else {
      attr(weights, "sum to one") <- TRUE
    }
  }

  return(list(id=id, weights=weights));
}

############################################################################
## HISTORY:
## 2015-08-03
## o Created.
############################################################################

