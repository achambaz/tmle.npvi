###########################################################################/**
# @RdocFunction validateArgumentObs
#
# @title "Validates argument 'obs'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations or indices (when 'family' is 'parsimonious').}
#   \item{allowIntegers}{A @logical (default 'TRUE'). If 'FALSE', 'obs[, "X"]' must have
#      positive mass in 0.}
# }
#
# \value{
#  Returns the (possibly slightly modified) input @matrix.
# }
#
# @author
#
#
#*/###########################################################################
validateArgumentObs <- function(obs, allowIntegers=TRUE) {
  ## if (!is.matrix(obs)) {
  ##   throw("Argument 'obs' should be a matrix");
  ## }
  if (!is.matrix(obs)) {
    throw("Argument 'obs' should be a matrix");
  }

  ## Ensuring that columns 'X' (exposure)  and 'Y' (outcome) are present
  varNames <- c("X", "Y");
  if (ncol(obs)==3) {
    ## Ensuring that column 'W' (covariate) is present too
    varNames <- c("W", varNames);
  }
  nms <- colnames(obs);
  m <- match(varNames, nms);
  if (nrow(obs)>0) {
    idxs <- which(is.na(m));
    if (length(idxs)) {
      throw("Missing column:", varNames[idxs]);
    }
    if (allowIntegers) {
      reject <- sum(obs[, "X"]==0)==0 & !is.integer(obs[, "X"])
    } else {
      reject <- is.integer(obs[, "X"]) | sum(obs[, "X"]==0)==0 
    }
    if (reject) {
      throw("Exposure (obs[,'X']==0) must have positive mass if 'tabulate' is FALSE");
      ## if 'tabulate' is TRUE, then this test fails to detect case
      ## where no index corresponding to a zero X appears in 'obs'
    }
  }
  obs;
}

############################################################################
## HISTORY:
## 2011-04-22
## o Created.
############################################################################

