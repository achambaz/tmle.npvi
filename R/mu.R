###########################################################################/**
# @RdocMethod getMumin
# @alias getMumin
#
# @title "Returns the value of parameter mumin"
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
#  Returns a @scalar, the value of parameter \var{mumin}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getMumin", "NPVI", function(this, ...) {
  this$.mumin;
})

###########################################################################/**
# @RdocMethod getMumax
# @alias getMumax
#
# @title "Returns the value of parameter mumax"
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
#  Returns a @scalar, the value of parameter \var{mumax}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getMumax", "NPVI", function(this, ...) {
  this$.mumax;
})


###########################################################################/**
# @RdocMethod getMu
# @alias getMu
#
# @title "Returns the value of parameter mu"
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
#  Returns a @function, the value of parameter \var{mu}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getMu", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.mu;
  } else {
    this$.mutab;
  } 
})



###########################################################################/**
# @RdocMethod setMu
# @alias setMu
#
# @title "Sets the value of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mu}{A @function that estimates the conditional expectation of
#       the DNA copy number given the DNA methylation level: E(X|W)}
# }
#
# \value{
#  Returns a NPVI @object containing the estimated \var{mu}, after
#  thresholding based on mumin and mumax. 
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getMu
# }
#
#*/###########################################################################
setMethodS3("setMu", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }
  
  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMu <- function(W) {
    threshold(mu(W), min=mumin, max=mumax)
  }
  this$.mu <- thresholdedMu ;
})



###########################################################################/**
# @RdocMethod setMuTab
# @alias setMuTab
#
# @title "Sets the value of the tabulated version of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mu}{A @function that estimates the conditional expectation of
#       the DNA copy number given the DNA methylation level: E(X|W)}
# }
#
# \value{
#  Returns a NPVI @object containing the tabulated version of the estimated \var{mu}, after
#  thresholding based on mumin and mumax.
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getMu
# }
#
#*/###########################################################################
setMethodS3("setMuTab", "NPVI", function(this, mu, ...) {
  ## Argument 'mu':
  if ((!is.null(mu))  && (mode(mu)!="function")) {
    throw("Argument \var{mu} should be of mode 'function', not ", mode(mu));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMu <- function(W) {
    threshold(mu(W), min=mumin, max=mumax)
  }

  this$.mutab <- thresholdedMu
})

###########################################################################/**
# @RdocMethod initializeMu
# @alias initializeMu
#
# @title "Initializes the estimation of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{muAux}{A @function, the estimated mu/(1-g).}
#   \item{g}{A @function, the estimated g:W -> P(X=0|W).}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing \var{mu=muAux*(1-g)}, after
#  thresholding based on mumin and mumax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateMu"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("initializeMu", "NPVI", function(this, muAux, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'muAux':
  if (mode(muAux) != "function") {
    throw("Argument 'muAux' should be a function, not a ", mode(muAux));
  }
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  mumin <- getMumin(this)
  mumax <- getMumax(this)
  thresholdedMuAux <- function(W) {
    threshold(muAux(W), min=mumin, max=mumax)
  }

  ## mu (after thresholding)
  mu <- function(W) {
    thresholdedMuAux(W)*(1-g(W))
  }
  setMu(this, mu)

  ## tabulated version of 'mu'
  fW <- getFW(this);
  obs <- getObs(this);
  MUTAB <- mu(fW(obs)); ## a *vector*, not a function
  mutab <- function(ii) {
    MUTAB[ii];
  }
  setMuTab(this, mutab)
})


###########################################################################/**
# @RdocMethod updateMu
# @alias updateMu
#
# @title "Updates the current estimation of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{mu}
#     at its current estimated value.}
#   \item{eps}{A @numeric, the "magnitude" of the increment to be performed.}
#   \item{exact}{If @TRUE, mu is updated using the exact relationship between
#     mu_epsilon and mu. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used.}
#   \item{effICW}{Conditional expectation of efficient influence curve given W,
#     only required when \var{exact} is @TRUE.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{mu}, after
#  thresholding based on mumin and mumax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateMu", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateMuNonTab(this, dev, exact=exact, effICW, ...)
  updateMuTab(this, dev, exact=exact, effICW, ...)
})

###########################################################################/**
# @RdocMethod updateMuNonTab
# @alias updateMuNonTab
#
# @title "Updates the non-tabulated version of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{mu}
#     at its current estimated value.}
#   \item{eps}{A @numeric, the "magnitude" of the increment to be performed.}
#   \item{exact}{If @TRUE, mu is updated using the exact relationship between
#     mu_epsilon and mu. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used.}
#   \item{effICW}{Conditional expectation of efficient influence curve given W,
#     only required when \var{exact} is @TRUE.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{mu}, after
#  thresholding based on mumin and mumax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateMuNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)

  mu <- getMu(this, tabulate=FALSE);
  
  if (!exact) { ## if do not use exact expression
    mu1 <- function(W) {
      mu(W) + eps * dev(W);
    }
  } else { ## if use exact expression
    mu1 <- function(W) {
      muW <- mu(W);
      theEffICW <- effICW(W)
      numerator <- muW + eps * (dev(W) + muW*theEffICW);
      denominator <- 1 + eps*theEffICW;
      numerator/denominator;
    }
  }
  setMu(this, mu1);
})

###########################################################################/**
# @RdocMethod updateMuTab
# @alias updateMuTab
#
# @title "Updates the tabulated version of parameter mu"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{mu}
#     at its current estimated value.}
#   \item{eps}{A @numeric, the "magnitude" of the increment to be performed.}
#   \item{exact}{If @TRUE, mu is updated using the exact relationship between
#     mu_epsilon and mu. If @FALSE, a first-order Taylor expression wrt epsilon
#     is used.}
#   \item{effICW}{Conditional expectation of efficient influence curve given W,
#     only required when \var{exact} is @TRUE.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{mu}, after
#  thresholding based on mumin and mumax.
# }
#
# @author
#
# \seealso{
#   @seemethod "estimateEpsilon"
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateMuTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  if (mode(dev) != "function") {
    throw("Argument 'dev' should be a function, not a ", mode(dev));
  }

  ## Argument 'exact':
  exact <- Arguments$getLogical(exact);
  
  ## Argument 'effICW':
  if (exact) {
    if (missing(effICW)) {
      throw("Argument 'effICW' is required when 'exact' is TRUE");
    }
    if (mode(effICW) != "function") {
      throw("Argument 'dev' should be a function, not a ", mode(dev));
    }
  }

  fW <- getFW(this)
  eps <- getEpsilon(this)

  mu <- getMu(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  muW <- mu(obs[, "W"])
  rm(mu, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  rm(obs)
  
  if (!exact) { ## do not use exact expression
    mu1W <- muW + eps * devW;
  } else { ## use exact expression
    theEffICW <- effICW(W)
    ## the above should use the tabulated or real versions of mu and theta0
    ## depending on tabulate, because effICW works on true values or
    ## indices depending on 'tabulate' (see how it is created in 'NPVI.update')
    numerator <- muW + eps * (devW + muW*theEffICW);
    denominator <- 1 + eps*theEffICW;
    mu1W <- numerator/denominator;
  }

  mu1tab <- function(ii) {
    mu1W[ii]
  }
  setMuTab(this, mu1tab)
})

############################################################################
## HISTORY:
## 2011-09-27
## o Truncation for 'mu' now performed on 'muAux' directly.
## 2011-05-03
## o Made two different methods for updating tabulated and non-tabulated
##   versions of 'mu'.
## 2011-05-02
## o Created.
############################################################################

