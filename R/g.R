###########################################################################/**
# @RdocMethod getGmin
# @alias getGmin
#
# @title "Returns the value of parameter gmin"
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
#  Returns a @scalar, the value of parameter \var{gmin}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getGmin", "NPVI", function(this, ...) {
  this$.gmin;
})

###########################################################################/**
# @RdocMethod getGmax
# @alias getGmax
#
# @title "Returns the value of parameter gmax"
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
#  Returns a @scalar, the value of parameter \var{gmax}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getGmax", "NPVI", function(this, ...) {
  this$.gmax;
})


###########################################################################/**
# @RdocMethod getG
# @alias getG
#
# @title "Returns the value of parameter g"
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
#  Returns a @function, the value of parameter \var{g}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getG", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    this$.g;
  } else {
    this$.gtab;
  }
})

###########################################################################/**
# @RdocMethod setG
# @alias setG
#
# @title "Sets the value of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{g}{A @function that estimates the conditional probability of
#       neutral copy number (X==0) given the DNA methylation level: P(X=0|W)}
# }
#
# \value{
#  Returns a NPVI @object containing the estimated \var{g}, after
#  thresholding based on gmin and gmax.
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getG
# }
#
#*/###########################################################################
setMethodS3("setG", "NPVI", function(this, g, ...) {
  ## Argument 'g':
  if ((!is.null(g))  && (mode(g)!="function")) {
    throw("Argument \var{g} should be of mode 'function', not ", mode(g));
  }

  gmin <- getGmin(this)
  gmax <- getGmax(this)
  thresholdedG <- function(W) {
    threshold(g(W), min=gmin, max=gmax)
  }
  this$.g <- thresholdedG
})

###########################################################################/**
# @RdocMethod setGTab
# @alias setGTab
#
# @title "Sets the value of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{g}{A @function that estimates the conditional probability of
#       neutral copy number (X==0) given the DNA methylation level: P(X=0|W)}
# }
#
# \value{
#  Returns a NPVI @object containing the tabulated version of the estimated \var{g}, after
#  thresholding based on gmin and gmax.
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getG
# }
#
#*/###########################################################################
setMethodS3("setGTab", "NPVI", function(this, g, ...) {
  ## Argument 'g':
  if ((!is.null(g))  && (mode(g)!="function")) {
    throw("Argument \var{g} should be of mode 'function', not ", mode(g));
  }

  gmin <- getGmin(this)
  gmax <- getGmax(this)
  thresholdedG <- function(W) {
    threshold(g(W), min=gmin, max=gmax)
  }

  this$.gtab <- thresholdedG
})


###########################################################################/**
# @RdocMethod initializeG
# @alias initializeG
#
# @title "Initializes the estimation of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{g}
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
#  Returns a NPVI @object containing the updated \var{g}, after
#  thresholding based on gmin and gmax.
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
setMethodS3("initializeG", "NPVI", function(this, g, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'g':
  if (mode(g) != "function") {
    throw("Argument 'g' should be a function, not a ", mode(g));
  }

  ## g
  setG(this, g)

  ## tabulated version of 'g'
  fW <- getFW(this);
  obs <- getObs(this);
  GTAB <- g(fW(obs)); ## a *vector*, not a function
  gtab <- function(ii) {
    GTAB[ii];
  }
  setGTab(this, gtab)
})

###########################################################################/**
# @RdocMethod updateG
# @alias updateG
#
# @title "Updates the current estimation of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{g}
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
#  Returns a NPVI @object containing the updated \var{g}, after
#  thresholding based on gmin and gmax.
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
setMethodS3("updateG", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
  updateGNonTab(this, dev, exact=exact, effICW, ...)
  updateGTab(this, dev, exact=exact, effICW, ...)
})

###########################################################################/**
# @RdocMethod updateGNonTab
# @alias updateGNonTab
#
# @title "Updates the non-tabulated version of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{g}
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
#  Returns a NPVI @object containing the updated \var{g}, after
#  thresholding based on gmin and gmax.
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
setMethodS3("updateGNonTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
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
  tabulate <- getTabulate(this)

  g <- getG(this, tabulate=FALSE)

  if (!exact) { ## do not use exact expression
    g1 <- function(W) {
      logit <- qlogis
      expit <- plogis
      
      gW <- g(W);
      res <- logit(gW) + eps * dev(W) * 1/(gW*(1-gW));
      expit(res);
    }
  } else { ## use exact expression
    g1 <- function(W) {
      gW <- g(W)
      theEffICW <- effICW(W)
      numerator <- gW + eps * (dev(W) + gW*theEffICW);
      denominator <- 1 + eps*theEffICW;
      out <- numerator/denominator;
      return(out)
    }
  }
  setG(this, g1);
})


###########################################################################/**
# @RdocMethod updateGTab
# @alias updateGTab
#
# @title "Updates the tabulated version of parameter g"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @function, the "derivative" of parameter \var{g}
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
#  Returns a NPVI @object containing the updated \var{g}, after
#  thresholding based on gmin and gmax.
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
setMethodS3("updateGTab", "NPVI", function(this, dev, exact=TRUE, effICW, ...) {
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
  tabulate <- getTabulate(this)

  g <- getG(this, tabulate=TRUE)
  obs <- getObs(this, tabulate=TRUE);
  gW <- g(obs[, "W"])
  rm(g, obs)

  obs <- getObs(this)
  devW <- dev(fW(obs))
  W <- obs[, "W"]
  rm(obs)
    
  if (!exact) { ## do not use exact expression
    logit <- qlogis
    expit <- plogis
      
    res <- logit(gW) + eps * devW * 1/(gW*(1-gW));
    g1W <- expit(res);
  } else { ## use exact expression
    theEffICW <- effICW(W)
    ## the above should use the tabulated or real versions of mu and theta0
    ## depending on tabulate, because effICW works on true values or
    ## indices depending on 'tabulate' (see how it is created in 'NPVI.update')
   
    numerator <- gW + eps * (devW + gW*theEffICW);
    denominator <- 1 + eps*theEffICW;
    g1W <- numerator/denominator;
  }

  g1tab <- function(ii) {
    g1W[ii]
  }
  setGTab(this, g1tab)
})

############################################################################
## HISTORY:
## 2011-05-03
## o Made two different methods for updating tabulated and non-tabulated
##   versions of 'g'.
## 2011-05-02
## o Created.
############################################################################

