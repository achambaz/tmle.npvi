###########################################################################/**
# @RdocFunction simulateData
#
# @title "Simulates observations (to be used in the estimation of \var{psi})"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{B}{A @numeric, number of observations to be generated}
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{g}{A @function that estimates the conditional probability of
#     neutral copy number (X==0) given the DNA methylation level: P(X=0|W)}
#   \item{mu}{A @function that estimates the conditional expectation of
#     the DNA copy number given the DNA methylation level: E(X|W)}
#   \item{sigma2}{A @numeric, the squared expectation of DNA copy number
#     relative to the neutral state: E(X^2).}
#   \item{family}{A @character, the distribution family of errors in the
#     simulation scheme. Either 'parsimonious' or  'gaussian'. If 'parsimonious' then
#     the conditional distribution of X given W and X!=0 is a mixture of Dirac masses
#     at three different observed Xs. If 'gaussian' then the latter conditional distribution
#     is Gaussian. In either cased, the conditional distributions meets the constraints
#     on the conditional mean and variance.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#   Returns a @matrix of observations with 3 columns:
#   \describe{
#    \item{Y}{expression level}
#    \item{X}{DNA copy number}
#    \item{W}{DNA methylation level.}
#   }
# }
#
# @author
#
# \seealso{
#   @see "getSample"
#   @see "NPVI"
# }
#
#*/###########################################################################
simulateData <- function(B, W, X, g, mu, sigma2, theta=NULL, Y=NA, weightsW=rep(1, length(W)), family=c("parsimonious", "gaussian"), verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'B':
  B <- Arguments$getNumeric(B);
  
  ## Argument 'W':
  W <- Arguments$getNumerics(W);

  ## Argument 'X':
  X <- Arguments$getNumerics(X);
   
  ## Argument 'g':
  mode <- mode(g);
  if (mode != "function") {
    throw("Argument 'g' should be of mode 'function', not '", mode);
  }

  ## Argument 'mu':
  mode <- mode(mu);
  if (mode != "function") {
    throw("Argument 'mu' should be of mode 'function', not '", mode);
  }

  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);
  if (is.na(sigma2)) {
    throw("Argument 'sigma2' should be provided...")
  }

  ## Argument 'theta'
  if (!is.null(theta)) {
    mode <- mode(theta)
    if (mode != "function") {
      throw("Argument 'theta' should be 'NULL' or of mode 'function', not '",  mode)
    }
  }
  
  ## Argument 'family':
  family <- match.arg(family);
  
  ## Argument 'Y'
  Y <- Arguments$getNumerics(Y);
  if (!is.null(theta)) {
    if (any(is.na(Y)) & family=="gaussian") {
      throw("Argument 'Y' of mode 'numerics' should be provided when 'family' is 'gaussian'.")
    }
  }
    
  ## Argument 'weightsW':
  weightsW <- Arguments$getNumerics(weightsW);
  nr <- length(W)
  if (length(weightsW) != nr) {
    throw("Length of 'weightsW' must match length of W!")
  }

  ## Argument 'verbose':
  verbose <- Arguments$getLogical(verbose);

  ## the probability P(X=0) directly computed on the whole dataset
  meanGW <- mean(X==0);
  whichXisZero <- which(X==0);
  whichXisNotZero <- which(X!=0);
  obsX <- X[whichXisNotZero];
  rm(X);

  WB <- sample(W, B, replace=TRUE, prob=weightsW)
  XB <- rep(NA, B)
  YB <- rep(NA, B)
  muWB <- mu(WB)
  gWB <- g(WB)
  ##
  U <- (runif(B) >= gWB)
  if (family=="gaussian") {
    XB[!U] <- 0
  } else if (family=="parsimonious") {
    XB[!U] <- whichXisZero[1] ## first index of row with X equal to 0
  }
  ##
  muW <- muWB[U]
  gW <- gWB[U]
  condMeanX <- muW/(1-gW)

  ##
  parameters <- list(meanGW=meanGW, muWB=muWB, gWB=gWB, U=U)
  if (family=="gaussian") {
    ## Note: here, 'tabulate' is necessarily FALSE
    ## hence 'W' are actual observations and not indices
    sigma2Bis <- sigma2 - mean(muWB^2/(1-gWB))
    if (sigma2Bis <=0) {
      cat("sigma2:\n")
      print(sigma2)
      cat("mean(muWB^2/(1-gWB)):\n")
      print(mean(muWB^2/(1-gWB)))
      throw("Problem here!... Parameter 'sigma2Bis' should be positive...")
    }
    condVarX <- sigma2Bis/(1-gW)
    XB[U] <- rnorm(sum(U), mean=condMeanX, sd=sqrt(condVarX))
    if (!is.null(theta)) {
      YB <- rnorm(B, mean=theta(cbind(X=XB, W=WB)), sd=sd(Y))
    }
  } else if (family=="parsimonious") {
    indices <- simulateParsimoniouslyXgivenW(WB[U], obsX, condMeanX, sigma2, parameters)
    XB[U] <- whichXisNotZero[indices]
    if (!is.null(theta)) {
      T <- theta(cbind(X=XB, W=WB))
      YB <- simulateParsimoniouslyYgivenXW(T, Y)
    }
  }
  obsB <- cbind(W=WB, X=XB, Y=YB)
  rm(WB, XB, YB, U)
  obsB
}

############################################################################
## HISTORY:
## 2011-04-08
## o Adapted for 'parsimonious' case;
##   in particular, changed arguments of function!
## 2011-03-23
## o Renamed Argument 'O' into 'obs'.
## 2010-12-08
## o Arguments are now validated.
## o Now Var(X|W) depends on W.  It is sigma2/1-g(0|W)
##   instead of sigma2/1-E[g(0|W)]
## 2010-11-26
## o Created.
############################################################################

