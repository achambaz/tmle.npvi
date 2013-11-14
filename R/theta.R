###########################################################################/**
# @RdocMethod getThetamin
# @alias getThetamin
#
# @title "Returns the value of parameter thetamin"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a @scalar, the value of parameter \var{thetamin}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getThetamin", "NPVI", function(this, ...) {
  this$.thetamin;
})

###########################################################################/**
# @RdocMethod getThetamax
# @alias getThetamax
#
# @title "Returns the value of parameter thetamax"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a @scalar, the value of parameter \var{thetamax}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getThetamax", "NPVI", function(this, ...) {
  this$.thetamax;
})


###########################################################################/**
# @RdocMethod getTheta
# @alias getTheta
#
# @title "Returns the value of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a @function, the value of parameter \var{theta}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getTheta", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.theta;
  } else {
    this$.thetatab;
  } 
})

###########################################################################/**
# @RdocMethod setTheta
# @alias setTheta
#
# @title "Sets the value of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{theta}{A @function that estimates the conditional expectation of
#    the expression level given DNA copy number and DNA methylation: E(Y|X,W).
# }
#
# \value{
#  Returns a NPVI @object containing the estimated \var{theta}, after
#  thresholding based on thetamin and thetamax. 
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getTheta
# }
#
#*/###########################################################################
setMethodS3("setTheta", "NPVI", function(this, theta, ...) {
  ## Argument 'theta':
  if ((!is.null(theta))  && (mode(theta)!="function")) {
    throw("Argument \var{theta} should be of mode 'function', not ", mode(theta));
  }
  
  thetamin <- getThetamin(this)
  thetamax <- getThetamax(this)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta <- function(XW) {
    threshold(theta(XW), min=thetamin, max=thetamax)
  }
  this$.theta <- thresholdedTheta;

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta0' accordingly
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta0 <- function(W) {
    XW <- cbind(X=0, W=W);
    threshold(theta(XW), min=thetamin, max=thetamax);
  }
  this$.theta0 <- thresholdedTheta0;
})




###########################################################################/**
# @RdocMethod setThetaTab
# @alias setThetaTab
#
# @title "Sets the value of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{theta}{A @function that estimates the conditional expectation of
#    the expression level given DNA copy number and DNA methylation: E(Y|X,W).}
# }
#
# \value{
#  Returns a NPVI @object containing the tabulated version of the estimated \var{theta}, after
#  thresholding based on thetamin and thetamax.
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getTheta
# }
#
#*/###########################################################################
setMethodS3("setThetaTab", "NPVI", function(this, theta, theta0, ...) {
  ## Argument 'theta':
  if ((!is.null(theta))  && (mode(theta)!="function")) {
    throw("Argument \var{theta} should be of mode 'function', not ", mode(theta));
  }

  ## Argument 'theta0':
  if ((!is.null(theta0))  && (mode(theta0)!="function")) {
    throw("Argument \var{theta0} should be of mode 'function', not ", mode(theta0));
  }
  
  thetamin <- getThetamin(this)
  thetamax <- getThetamax(this)

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta <- function(XW) {
    threshold(theta(XW), min=thetamin, max=thetamax)
  }
  this$.thetatab <- thresholdedTheta

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## set 'theta0'
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  thresholdedTheta0 <- function(W) {
    threshold(theta0(W), min=thetamin, max=thetamax)
  }
  this$.theta0tab <- thresholdedTheta0
})


###########################################################################/**
# @RdocMethod getTheta0
# @alias getTheta0
#
# @title "Returns the value of parameter theta0"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a @function, the value of parameter \var{theta0}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getTheta0", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.theta0;
  } else {
    this$.theta0tab;
  }
})

###########################################################################/**
# @RdocMethod initializeTheta
# @alias initializeTheta
#
# @title "Initializes the estimation of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{theta}{A @function that estimates the conditional expectation of
#    the expression level given DNA copy number and DNA methylation: E(Y|X,W).}
# }
#
# \value{
#  Returns a NPVI @object containing the input function \var{theta}, after
#  thresholding based on thetamin and thetamax.
# }
#
# \details{
#   Argument 'theta' is a function of the _true_ observations (X,W), not
#   a function of indices.
# }
#
# @author
#
# \seealso{
#   @seemethod "init"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("initializeTheta", "NPVI", function(this, theta, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'theta':
  if (mode(theta) != "function") {
    throw("Argument 'theta' should be a function, not a ", mode(theta));
  }

  ## theta
  setTheta(this, theta)

  ## tabulated version of 'theta'
  fW <- getFW(this);
  fX <- getFX(this);
  obs <- getObs(this);

  eg <- expand.grid(X=obs[, "X"], W=obs[, "W"])
  OBSTAB <- cbind(X=fX(eg), fW(eg))

  THETATAB <- theta(OBSTAB)
  THETATAB <- matrix(THETATAB, nrow=nrow(obs), ncol=nrow(obs))
  thetatab <- function(xiwj) {
    stopifnot(is.matrix(xiwj) && ncol(xiwj)==2 && is.integer(xiwj))
    THETATAB[xiwj]
  }
  ## tabulated version of 'theta0'
  THETA0TAB <- theta(cbind(X=0, W=fW(obs)))  ## a vector, not a matrix !
  theta0tab <- function(wi) {
    ## stopifnot(is.matrix(wi) && ncol(wi)==1 && is.integer(wi))
    THETA0TAB[wi]
  }
  setThetaTab(this, thetatab, theta0tab);
})



###########################################################################/**
# @RdocMethod updateTheta
# @alias updateTheta
#
# @title "Updates the current value of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{theta}
#     at its current estimated value.}
#   \item{cleverCovTheta}{If @TRUE, theta is updated using the corresponding
#     "clever covariate".}
#   \item{exact}{If @TRUE, theta is updated using the exact relationship between
#     theta_epsilon and theta. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used. If \var{cleverCovTheta} is @TRUE, then by definition \var{exact} is
#     forced to @TRUE in the updating process of \var{theta}.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{theta}, after
#  thresholding based on thetamin and thetamax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "updateEstimation"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateTheta", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  updateThetaNonTab(this, dev, cleverCovTheta, exact=exact, ...)
  updateThetaTab(this, dev, cleverCovTheta, exact=exact, ...)
})

###########################################################################/**
# @RdocMethod updateThetaNonTab
# @alias updateThetaNonTab
#
# @title "Updates the non-tabulated version of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{theta}
#     at its current estimated value.}
#   \item{cleverCovTheta}{If @TRUE, theta is updated using the corresponding
#     "clever covariate".}
#   \item{exact}{If @TRUE, theta is updated using the exact relationship between
#     theta_epsilon and theta. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used. If \var{cleverCovTheta} is @TRUE, then by definition \var{exact} is
#     forced to @TRUE in the updating process of \var{theta}.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{theta}, after
#  thresholding based on thetamin and thetamax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "updateEstimation"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateThetaNonTab", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if ((mode(dev) != "function") & (!cleverCovTheta)) {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  ## Argument 'cleverCovTheta':
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);

  fW <- getFW(this)
  fX <- getFX(this)
  devThetaAux <- dev

  theta <- getTheta(this, tabulate=FALSE);
  
  if (!cleverCovTheta) {
    eps <- getEpsilon(this);
    g <- getG(this, tabulate=FALSE);
    mu <- getMu(this, tabulate=FALSE);
    sigma2 <- getSigma2(this);
    psi <- getPsi(this);
    if (!exact) { ## do not use clever covariate nor exact expression
      dev <- function(XW) {
        X <- XW[, 1, drop=TRUE];
        W <- XW[, 2, drop=FALSE];
        devThetaAux(XW) * (X - mu(W)/g(W)*(X==0))/sigma2;
      }
      theta1 <- function(XW) {
        theta(XW) + eps * dev(XW);
      }
    } else { ## do not use clever covariate, but exact expression
      theta1 <- function(XW) {
        X <- XW[, 1, drop=TRUE];
        W <- XW[, 2, drop=FALSE];
        TXW <- theta(XW);
        term1 <- devThetaAux(XW) * (X - mu(W)/g(W)*(X==0))/sigma2;
        term2 <- (X*(TXW-theta(cbind(X=0, W=W))-X*psi))/sigma2;
        numerator <- TXW + eps*(term1+term2*TXW);
        denominator <- 1 + eps*term2;
        numerator/denominator;
      }
    }
  } else { ## use clever covariate (hence, exact expression)
    if (!exact) {
      throw("Can't use a clever covariate for updating 'theta' if 'exact==FALSE' !");
    }
    eps <- getEpsilonTheta(this);
    dev <- getHTheta(this, tabulate=FALSE);
    theta1 <- function(XW) {
      theta(XW) + eps * dev(XW);
    }
  }
  setTheta(this, theta1)
})


###########################################################################/**
# @RdocMethod updateThetaTab
# @alias updateThetaTab
#
# @title "Updates the tabulated version of parameter theta"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{theta}
#     at its current estimated value.}
#   \item{cleverCovTheta}{If @TRUE, theta is updated using the corresponding
#     "clever covariate".}
#   \item{exact}{If @TRUE, theta is updated using the exact relationship between
#     theta_epsilon and theta. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used. If \var{cleverCovTheta} is @TRUE, then by definition \var{exact} is
#     forced to @TRUE in the updating process of \var{theta}.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{theta}, after
#  thresholding based on thetamin and thetamax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "updateEstimation"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateThetaTab", "NPVI", function(this, dev, cleverCovTheta, exact=TRUE, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if ((mode(dev) != "function") & (!cleverCovTheta)) {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  
  ## Argument 'cleverCovTheta':
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);

  fW <- getFW(this)
  fX <- getFX(this)
  devThetaAux <- dev

  theta <- getTheta(this, tabulate=TRUE);
  theta0 <- getTheta0(this, tabulate=TRUE)

  obs <- getObs(this, tabulate=TRUE);
  nr <- nrow(obs)

  ## rm(obs)
  ## XW <- as.matrix(expand.grid(X=1:nr, W=1:nr))
  ## thetaXW <- matrix(theta(XW), nrow=nr, ncol=nr)

  XW <- as.matrix(expand.grid(X=obs[, "X"], W=obs[, "W"]))
  thetaXW <- matrix(theta(XW), nrow=nr, ncol=nr)
  rm(obs)

  if (!cleverCovTheta) {
    eps <- getEpsilon(this);
    g <- getG(this, tabulate=TRUE);
    mu <- getMu(this, tabulate=TRUE);
    sigma2 <- getSigma2(this);
    psi <- getPsi(this);

    obs <- getObs(this, tabulate=TRUE);
    gW <- g(obs[, "W"])
    gW <- t(matrix(gW, nrow=nr, ncol=nr))
    muW <- mu(obs[, "W"])
    muW <- t(matrix(muW, nrow=nr, ncol=nr))

    rm(g, mu, obs)

    obs <- getObs(this)
    eg <- expand.grid(X=obs[, "X"], W=obs[, "W"])
    OBSTAB <- cbind(X=fX(eg), fW(eg))
    X <- matrix(fX(obs), nrow=nr, ncol=nr)
    rm(obs)
    
    devThetaAuxXW <- matrix(devThetaAux(OBSTAB), nrow=nr, ncol=nr)
    devXW <- devThetaAuxXW * (X - muW/gW*(X==0))/sigma2;
    
    if (!exact) { ## do not use clever covariate nor exact expression
      theta1XW <- thetaXW + eps*devXW
    } else { ## do not use clever covariate, but exact expression
      theta0W <- t(matrix(theta0(1:nr), nrow=nr, ncol=nr))
      term2 <- ( X*(thetaXW - theta0W -X*psi) )/sigma2;
      numerator <- thetaXW + eps*(devXW+term2*thetaXW);
      denominator <- 1 + eps*term2;
      theta1XW <- numerator/denominator;
    }
  } else { ## use clever covariate (hence, exact expression)
    if (!exact) {
      throw("Can't use a clever covariate for updating 'theta' if 'exact==FALSE' !");
    }
    eps <- getEpsilonTheta(this);
    dev <- getHTheta(this, tabulate=TRUE);
    devXW <- dev(expand.grid(X=1:nr, W=1:nr))
    devXW <- matrix(devXW, nrow=nr, ncol=nr)
    theta1XW <- thetaXW + eps*dev(XW)
  }
  ## retrieving 'theta01W' from 'theta1XW'
  obs <- getObs(this, tabulate=FALSE)
  idx <- which(obs[, "X"]==0)[1]
  theta01W <- theta1XW[idx, ]

  ## formatting
  THETA1TAB <- theta1XW
  theta1tab <- function(xiwj) {
    stopifnot(is.matrix(xiwj) && ncol(xiwj)==2 && is.integer(xiwj))
    THETA1TAB[xiwj]
  }
  THETA01TAB <- theta01W
  theta01tab <- function(wi) {
    ## stopifnot(is.matrix(wi) && ncol(wi)==1 && is.integer(wi))
    THETA01TAB[wi]
  }
  
  setThetaTab(this, theta1tab, theta01tab);
})


############################################################################
## HISTORY:
## 2011-05-03
## o Made two different methods for updating tabulated and non-tabulated
##   versions of 'theta'.
## 2011-05-02
## o Created.
############################################################################

