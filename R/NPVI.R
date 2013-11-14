###########################################################################/**
# @RdocClass NPVI
#
# @title "The NPVI class"
#
# \description{
#  @classhierarchy
#
# This class represents parameter estimates of the relationship between 
# copy number (X) and expression (Y), accounting for DNA methylation (W)
# in terms of non-parametric variable importance (NPVI).
#
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{family}{A @character, describing how X should be conditionally simulated
#     given W. Currently two options: 'parsimonious' and 'gaussian'.}
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \details{
#   Objects of class NPVI contain the following slots:
#   \describe{
#     \item{g}{A @function that estimates the conditional probability of
#       neutral copy number (X==0) given the DNA methylation level: P(X=0|W)}
#     \item{mu}{A @function that estimates the conditional expectation of
#       the DNA copy number given the DNA methylation level: E(X|W)}
#     \item{theta}{A @function that estimates the conditional expectation of
#       the expression level given DNA copy number and DNA methylation: E(Y|X,W)}
#     \item{sigma2}{A @numeric, the squared expectation of DNA copy number
#       relative to the neutral state: E(X^2).}
#     \item{psi}{A @numeric, the parameter estimate}
#     \item{psiPn}{A @numeric, an alternative parameter estimate (invalid until final step).}
#     \item{gmin, gmax, mumin, mumax, thetamin, thetamax}{Thresholds to pass to 'threshold'.}
#     \item{efficientInfluenceCurve}{The estimated efficient influence
#       curve of parameter \var{psi}, in the form of a @matrix with 3
#       columns: the two components of the curve, and their sum.}
#     \item{epsilon}{A @numeric @vector of length 2, the estimated magnitude
#       of parameter updates that has been used for \acronym{TML} estimation.}
# }
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/########################################################################### 

setConstructorS3("NPVI", function(obs=matrix(nrow=0, ncol=3, dimnames=list(NULL, c("W", "X", "Y"))),
                                  f=identity,
                                  gmin=0.01, gmax=1-gmin, mumin=-Inf,
                                  mumax=Inf, thetamin=-Inf, thetamax=Inf,
                                  family=c("parsimonious", "gaussian"), tabulate=TRUE,
                                  stoppingCriteria=list(mic=0.03, div=0.01, psi=0.1),
                                  conf.level=0.95,
                                  ...,
                                  verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- validateArgumentObs(obs, allowIntegers=FALSE);
  
  ## Argument 'f':
  if (!((mode(f)=="function") && (f(0)==0))) {
    throw("Argument 'f' must be a function such that f(0)=0.")
  }

  ## Arguments 'gmin' and 'gmax':
  gmin <- Arguments$getNumeric(gmin);
  gmax <- Arguments$getNumeric(gmax);
  if (gmin>=1 | gmin<=0) {
    throw("Argument 'gmin' must belong to ]0,1[");
  }
  if (gmax>=1 | gmax<=0) {
    throw("Argument 'gmax' must belong to ]0,1[");
  }
  if (gmin>=gmax) {
    throw("Argument 'gmin' must be smaller than argument 'gmax'.")
  }

  ## Arguments 'mumin' and 'mumax':
  mumin <- Arguments$getNumeric(mumin);
  mumax <- Arguments$getNumeric(mumax);
  if (mumin>=mumax) {
    throw("Argument 'mumin' must be smaller than argument 'mumax'.")
  }

  ## Arguments 'thetamin' and 'thetamax':
  thetamin <- Arguments$getNumeric(thetamin);
  thetamax <- Arguments$getNumeric(thetamax);
  if (thetamin>=thetamax) {
    throw("Argument 'thetamin' must be smaller than argument 'thetamax'.")
  }

  ## Argument 'family'
  family <- match.arg(family)

  ## Argument 'tabulate'
  tabulate <- Arguments$getLogical(tabulate)

  ## Argument 'stoppingCriteria'
  mic.tol <- Arguments$getNumeric(stoppingCriteria$mic)
  attr(mic.tol, "label") <- "scaled empirical mean of estimating function"
  div.tol <- Arguments$getNumeric(stoppingCriteria$div)
  attr(div.tol, "label") <- "TV distance between P_n^k and P_n^{k+1}"
  psi.tol <- Arguments$getNumeric(stoppingCriteria$psi)
  attr(psi.tol, "label") <- "change between successive values of \"psi\""

  ## Argument 'conf.level'
  conf.level <- Arguments$getNumeric(conf.level, c(0, 1))

  ## Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  } 

  if (tabulate & family=="gaussian") {
    throw("Unauthorized, because cannot tabulate all functions if family is 'gaussian'.")
  }
  if (!tabulate & family=="parsimonious") {
    throw("Unauthorized, because not implemented (it is a sub-optimal combination of options).")
  }

  obs[, "X"] <- f(obs[, "X"])
  
  theW <- setdiff(colnames(obs), c("X", "Y"))
  if (!tabulate & length(theW)>1) {
    throw("Multivariate 'W' handled only if 'tabulate' is TRUE")
  }
  
  nms <-  c("eps", "lli", "mic1", "epsT", "lliT", "mic2", "psi", "psi.sd", "psiPn", "psiPn.sd", "mic", "div", "sic", "phi", "sicAlt")
  history <- matrix(NA, 0, length(nms));
  colnames(history) <- nms

  conv <- NA
  attr(conv, "msg") <- character(0)
  
  extend(Object(), "NPVI",
         .obs=obs, #.flavor=flavor,
         .g=NULL, .mu=NULL, .theta=NULL, .theta0=NULL, .weightsW=rep(1, nrow(obs)), 
         .gtab=NULL, .mutab=NULL, .thetatab=NULL, .theta0tab=NULL,
         .sigma2=NA, .psi=NA, .psi.sd=NA, .psiPn=NA, .psiPn.sd=NA,
         .gmin=gmin, .gmax=gmax, .mumin=mumin, .mumax=mumax,
         .thetamin=thetamin, .thetamax=thetamax,
         .family=family, .tabulate=tabulate, .epsilon=NA,
         .epsilonTheta=NA, .logLikIncr=NA, .logLikIncrTheta=NA, .div=NA,
         .efficientInfluenceCurve=matrix(NA, 0, 3), .history=history, .step=0,
         .stoppingCriteria=list(mic=mic.tol, div=div.tol, psi=psi.tol),
         .conv=conv, .conf.level=conf.level
         );
})


###########################################################################/**
# @RdocMethod getFlavor
# @alias getFlavor
#
# @title "Returns the estimation 'flavor' of the NPVI @object"
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
#  Returns a @character @value, the 'flavor' used for the esimation.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getFlavor", "NPVI", function(this, ...) {
  this$.flavor;
})

setMethodS3("getMicTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$mic;
})

setMethodS3("getDivTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$div;
})

setMethodS3("getPsiTol", "NPVI", function(this, ...) {
  this$.stoppingCriteria$psi;
})

setMethodS3("getConfLevel", "NPVI", function(this, ...) {
  this$.conf.level;
})

setMethodS3("getHistory", "NPVI", function(
  ### Returns the 'history' of the TMLE procedure.
  this,
### An object of class \code{TMLE.NPVI}.
  ...
### Not used.
  ) {
  ##alias<< getHistory
  ##seealso<< tmle.npvi
  this$.history;
  ###  Returns a \code{numeric}  \code{matrix} which encapsulates a summary of
  ###     the  TMLE procedure.  If \eqn{k} successive  updates were performed,
  ###       then   the   \code{matrix}   has   either   \eqn{k+1}   rows   (if
  ###      \code{cleverCovTheta}  was  set  to  \code{FALSE} in  the  call  to
  ###       \code{tmle.npvi})   or    \code{2k+1}   rows   (otherwise).    The
  ###     \code{matrix} has 14 columns:
  ### \itemize{
  ###       \item{\code{"eps"},   values of  the unique
  ###        fluctuation  parameter   (if  \code{cleverCovTheta}  was  set  to
  ###        \code{FALSE} in  the call  to \code{tmle.npvi}),  or  
  ###        values of the  parameter involved in the fluctuation of
  ###        the   joint  distribution  of  \eqn{(X,W)}   during  each  update
  ###       (otherwise).  }
  ###     \item{\code{"lli"},     increases in likelihood
  ###      yielded  by  each  update  (if  \code{cleverCovTheta}  was  set  to
  ###      \code{FALSE}  in the  call  to  \code{tmle.npvi}),  or 
  ###      increases in likelihood yielded by the fluctuation of the
  ###     joint distribution of \eqn{(X,W)} during each update (otherwise).}
  ###        \item{\code{"mic1"},   empirical means  of the  first
  ###     component of the efficient influence  curve at each step of the TMLE
  ###     procedure.}
  ### \item{\code{"epsT"},   values of  the fluctuation
  ### parameter involved in the fluctuation of the conditional distribution of
  ### \eqn{Y}  given \eqn{(X,W)} during each  update (if \code{cleverCovTheta}
  ### was  set to \code{TRUE} in  the call to  \code{tmle.npvi}), or \code{NA}
  ### (otherwise).}
  ###     \item{\code{"lliT"},  successive increases in likelihood
  ###       yielded by  the  fluctuation of  the  conditional distribution  of
  ###         \eqn{Y}    given   \eqn{(X,W)}    during    each   update    (if
  ###       \code{cleverCovTheta}  was  set  to  \code{TRUE} in  the  call  to
  ###      \code{tmle.npvi}), or \code{NA} (otherwise).}
  ### \item{\code{"mic2"},  empirical means  of the  second
  ###     component of the efficient influence  curve at each step of the TMLE
  ###     procedure.}
  ###   \item{\code{"psi"},   increasingly targeted
  ###     estimators \eqn{\Psi(P_n^k)} of  the  parameter  of  interest. The last one is the TMLE.  Their  computation  
  ###       involves  simulation  of  \code{B}  iid  copies  of
  ###    \eqn{(X,W)} under \eqn{P_n^k}. }
  ### \item{\code{"psi.sd"},  estimated standard deviations of
  ###  the   increasingly targeted  estimators of  the  parameter of
  ###  interest. The last one corresponds to the TMLE. The  computation involves  the  same \code{B}  iid copies  of
  ### \eqn{(X,W)} as above.}
  ###  \item{\code{"psiPn"},  same as \code{"psi"} except  that the *observed*
  ###   \eqn{(X_i,W_i)}  are  used  instead  of simulated  copies  drawn  from
  ###  \eqn{P_n^k}. Of course, \code{"psi"} must be favored.}
  ###   \item{\code{"psiPn.sd"},  same  as  \code{"psi.sd"}  except  that  the
  ###  *observed*  \eqn{(X_i,W_i)} are used instead of  simulated copies drawn
  ###  from \eqn{P_n^k}. Of course, \code{"psi.sd"} must be favored.}
  ###  \item{\code{"mic"},   empirical  means  of  the  efficient
  ### influence curve at each step of the TMLE procedure. This column is the sum of the \code{"mic1"} and \code{"mic2"} columns.}
  ### \item{\code{"div"},  total variation  distances between
  ### each pair  of successive distributions constructed in  the course of the
  ### TMLE procedure. }
  ### \item{\code{"sic"},  estimated standard deviations   of  the  efficient
  ### influence curve at each step of the TMLE procedure.}
  ###\item{\code{"phi"},    non-parametric     substitution    estimator    of
  ###\eqn{\phi=\Phi(P)}             where            \deqn{\Phi(P)            =
  ###\frac{E_P[f(X)Y]}{E_P[f(X)^2]},}{\Phi(P)  =  E_P[f(X)Y]  /  E_P[f(X)^2],}
  ###with  \eqn{P}  the  distribution   of  the  random  vector  \eqn{(W,X,Y)}. The alternative parameter \eqn{\phi} should be interpreted as the counterpart of \eqn{\psi} which neglects \eqn{W}. }
  ### \item{\code{"sicAlt"},  estimated standard deviations   of  the  efficient
  ### influence curve of \eqn{\Psi - \Phi} at each step of the TMLE procedure.}
###   }
})

###########################################################################/**
# @RdocMethod updateHistory
# @alias updateHistory
#
# @title "Updates the 'history' of the NPVI @object"
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
#  An object of class NPVI where one column containing the current values of
#  \var{epsilon}, \var{logLikIncr}, \var{epsilonTheta},
#  \var{logLikIncrTheta}, \var{psi} has been added to the 'history' @matrix.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateHistory", "NPVI", function(this, ...) {
  history <- getHistory(this);
  psi <- getPsi(this);
  psi.sd <- getPsiSd(this);  
  epsilon <- getEpsilon(this);
  logLikIncr <- getLogLikIncr(this);
  epsilonTheta <- getEpsilonTheta(this);
  logLikIncrTheta <- getLogLikIncrTheta(this);
  eic <- getEfficientInfluenceCurve(this);
  mic <- apply(eic, 2, mean);
  sic <- getSic(this);
  div <- getDivergence(this)
  psiPn <- getPsiPn(this);
  psiPn.sd <- getPsiPnSd(this);

  phi <- getPhi(this)
  sicAlt <- getSicAlt(this)
  
  currStep <- c(epsilon, logLikIncr, mic[1], epsilonTheta, logLikIncrTheta, mic[2], psi, psi.sd, psiPn, psiPn.sd, mic[3], div, sic, phi, sicAlt);
  step <- getStep(this);
  rownames.history <- rownames(history);
  history <- rbind(history, matrix(currStep, 1, length(currStep)));
  rownames(history) <- c(rownames.history,
                         paste("step", as.character(step), sep=""));
  
  this$.history <- history;
})

setMethodS3("getPsi", "NPVI", function(
### Returns the current value of the estimator.
    this,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< getPsi
  ##seealso<< tmle.npvi, getHistory, getPsiSd
  this$.psi;
  ### Retrieves  the current value  of the estimator \eqn{\Psi(P_n^k)}  of the
  ### parameter  of interest. Its computation involves  simulation of a large number of 
  ### iid copies of \eqn{(X,W)} under \eqn{P_n^k}.
})

setMethodS3("getPsiSd", "NPVI", function(
### Returns the current value of the estimated standard deviation of the current estimator. 
    this,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< getPsiSd
  ##seealso<< tmle.npvi, getHistory, getPsi
  this$.psi.sd;
  ### Retrieves  the estimated standard deviation of the current estimator \eqn{\Psi(P_n^k)}  of the
  ### parameter  of interest. Its computation involves  simulation of a large number of 
  ### iid copies of \eqn{(X,W)} under \eqn{P_n^k}.  
})

setMethodS3("getPhi", "NPVI", function(this, ...) {
  obs <- getObs(this)
  fX <- getFX(this)
  fY <- getFY(this)
  X <- fX(obs)
  Y <- fY(obs)
  sX2 <- mean(X^2)
  ## phi: estimator of phi_0 
  mean(X*Y)/sX2
})

setMethodS3("getSic", "NPVI", function(this, ...) {
  eic <- getEfficientInfluenceCurve(this);
  mic <- apply(eic, 2, mean);
  sd(eic[, 3])
})

setMethodS3("getSicAlt", "NPVI", function(this, ...) {
  obs <- getObs(this)
  fX <- getFX(this)
  fY <- getFY(this)
  X <- fX(obs)
  Y <- fY(obs)
  eic <- getEfficientInfluenceCurve(this)
  ## phi: estimator of phi_0 
  phi <- getPhi(this)
  ## sicAlt: estimated standard deviation to perform test of "psi_0 = phi_0"
  sX2 <- mean(X^2)
  infCurvePhi <- (X*Y-phi*X^2)/sX2
  sicAlt <- sd(eic[, 3]-infCurvePhi)
  sicAlt
})


###########################################################################/**
# @RdocMethod getPsiPn
# @alias getPsiPn
#
# @title "Returns an alternative estimated value of psi"
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
#  Returns a @numeric, the value of the 'estimate' of \var{psi} such that
#  the empirical measure applied to the first component of the efficient
#  influence curve at the current P equal zero. 
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getPsiPn", "NPVI", function(this, name, ...) {
  this$.psiPn;
})

setMethodS3("getPsiPnSd", "NPVI", function(this, name, ...) {
  this$.psiPn.sd;
})

###########################################################################/**
# @RdocMethod getStep
# @alias getStep
#
# @title "Returns the 'step' of the TMLE NPVI procedure."
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
#  Returns a @integer, the 'step' of the TMLE NPVI procedure.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getStep", "NPVI", function(this, name, ...) {
  this$.step;
})

###########################################################################/**
# @RdocMethod getDivergence
# @alias getDivergence
#
# @title "Returns  the 'divergence'  of the current  estimated data-generating
#  distribution   relative   to   the   previous   estimated   data-generating
#  distribution, expressed in terms of total variation.
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
#  Returns a @numeric, the 'divergence'.  
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getDivergence", "NPVI", function(this, name, ...) {
  this$.div;
})

###########################################################################/**
# @RdocMethod getDivergence
# @alias getDivergence
#
# @title "Returns  the 'divergence'  of the current  estimated data-generating
#  distribution   relative   to   the   previous   estimated   data-generating
#  distribution, expressed in terms of total variation.
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
#  Returns a @numeric, the 'divergence'.  
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("setDivergence", "NPVI", function(this, div, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'div':
  div <- Arguments$getNumeric(div)
  
  this$.div <- div;
})

setMethodS3("setConfLevel", "NPVI", function(
### Sets the confidence level of a \code{TMLE.NPVI} object.
    this,
### An object of class \code{TMLE.NPVI}.
    confLevel,
### A \code{numeric}, confidence interval level.
    ...
### Not used.
    ) {
  ##alias<< setConfLevel
  ##seealso<< as.character.NPVI
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'confLevel':
  conf.level <- Arguments$getNumeric(confLevel, range=c(0, 1))
  
  this$.conf.level <- conf.level;
})


###########################################################################/**
# @RdocMethod getFW
# @alias getFW
#
# @title "Builds wrapper function 'fW'"
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
#  Returns a @function, which returns true observations 'W'.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getFW", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    ## 'W' is necessarily unidimensional
    fW <- function(obs) {
      obs[, "W", drop=FALSE];
    }
  } else {
    ## 'W' needs not be unidimensional
    obsT <- getObs(this, tabulate=FALSE);
    fW <- function(obs) {
      ## obsT[obs[, "W"], "W", drop=FALSE];
      extractW(obsT[obs[, "W"], ])
    }
  }
  fW
})


###########################################################################/**
# @RdocMethod getFX
# @alias getFX
#
# @title "Builds wrapper function 'fX'"
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
#  Returns a @function, which returns true observations 'X'.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getFX", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    fX <- function(obs) {
      obs[, "X"];
    }
  } else {
    obsT <- getObs(this, tabulate=FALSE);
    fX <- function(obs) {
      obsT[obs[, "X"], "X"];
    }
  }
  fX
})

###########################################################################/**
# @RdocMethod getFY
# @alias getFY
#
# @title "Builds wrapper function 'fY'"
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
#  Returns a @function, which returns true observations 'Y'.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getFY", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  if (!tabulate) {
    fY <- function(obs) {
      obs[, "Y"];
    }
  } else {
    obsT <- getObs(this, tabulate=FALSE);
    fY <- function(obs) {
      obsT[obs[, "Y"], "Y"];
    }
  }
  fY
})




setMethodS3("getObs", "NPVI", function(
### Retrieves the \code{matrix} of observations involved in the TMLE procedure.
    this,
### An object of class \code{TMLE.NPVI}.
    tabulate,
### A \code{logical}, to  specify whether it is the original  data set that is
### retrieved (if \code{FALSE}) or a  tabulated version of it (otherwise), for
### internal use only.  If \code{tabulate} is missing then  the value attached
### to the input object is used.
    ...
### Not used.
    ) {
  ##alias<< getObs
  ##seealso<< tmle.npvi
  if (missing(tabulate)) {
    tabulate <- getTabulate(this);
  }
  if (!tabulate) {
    obs <- this$.obs;
  } else {
    n <- nrow(this$.obs)
    obs <- matrix(1:n, nrow=n, ncol=3);
    colnames(obs) <- c("W", "X", "Y");
  }
  obs
### Either the original data set involved in the TMLE procedure or a tabulated
### version of it.
})






###########################################################################/**
# @RdocMethod getFamily
# @alias getFamily
#
# @title "Returns the 'family' of the procedure"
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
#  Returns a @character, the value of parameter \var{family}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getFamily", "NPVI", function(this, ...) {
  this$.family;
})

###########################################################################/**
# @RdocMethod setFamily
# @alias setFamily
#
# @title "Sets the value of the 'family' parameter of the procedure"
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
#  Set the @character value of parameter \var{family}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("setFamily", "NPVI", function(this, family=c("parsimonious", "gaussian"), ...) {
  ## Argument
  family <- match.arg(family);
  this$.family <- family;
})


###########################################################################/**
# @RdocMethod getTabulate
# @alias getTabulate
#
# @title "Returns the 'tabulate' parameter."
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
#  Returns a @logical, indicating whether resort to tabulating or not.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getTabulate", "NPVI", function(this, ...) {
  this$.tabulate;
})




###########################################################################/**
# @RdocMethod getHTheta
# @alias getHTheta
#
# @title "Returns the function HTheta"
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
#  Returns a @function, the value of parameter \var{HTheta}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getHTheta", "NPVI", function(this, tabulate, ...) {
  if (missing(tabulate)) {
    tabulate <- getTabulate(this)
  }
  mu <- getMu(this, tabulate)
  g <- getG(this, tabulate)
  sigma2 <- getSigma2(this)
  fX <- getFX(this, tabulate);
    
  HTheta <- function(XW) {
    X <- fX(cbind(XW, Y=NA))
    W <- XW[, "W"]
    (X - (X==0)*mu(W)/g(W))/sigma2
  }
  HTheta;
})

###########################################################################/**
# @RdocMethod getSigma2
# @alias getSigma2
#
# @title "Returns the value of parameter sigma2"
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
#  Returns a @numeric, the value of parameter \var{sigma2}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getSigma2", "NPVI", function(this, ...) {
  this$.sigma2;
})

###########################################################################/**
# @RdocMethod setSigma2
# @alias setSigma2
#
# @title "Sets the value of parameter sigma2"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{sigma2}{A @numeric, the squared expectation of DNA copy number
#     relative to the neutral state: E(X^2).}
# }
#
# \value{
#  Returns a NPVI @object containing the estimated \var{sigma2}.
# }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getSigma2
# }
#
#*/###########################################################################
setMethodS3("setSigma2", "NPVI", function(this, sigma2, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2)

  this$.sigma2 <- sigma2;
})

setMethodS3("as.character", "NPVI", function(
### Returns a short string describing the NPVI object.
    x,
### An object of class \code{TMLE.NPVI}.
    ...
### Not used.
    ) {
  ##alias<< as.character
  ##seealso<< tmle.npvi
  this <- x;  ## To please R CMD check
  s <- sprintf("%s object:", class(this)[1]);
  s <- c(s, "")

  ## sample size
  s <- c(s, sprintf("Sample size: %s", nrow(getObs(this))))
  s <- c(s, "")
  
  ## psi
  s <- c(s, sprintf("Estimator of psi:\t\t%s", signif(getPsi(this), 3)));

  ## std error
  sic <- getSic(this)
  s <- c(s, sprintf("Estimated standard error:\t%s", signif(sic, 3)))
  s <- c(s, "")
  
  ## number of iterations
  step <- getStep(this)
  
  ## convergence ?
  conv <- getConv(this)
  if (!is.na(conv)) {
    if (conv) {
      msg <- paste("Convergence reached after", step, "iteration(s) because:")
      msg <- paste("\n", msg, "\n\t",
                   paste(attr(conv, "msg"), collapse="\n\t"),
                   sep="")
    } else {
      msg <- paste("Convergence not reached after", step, "iteration(s)")
    }

    mic <- getMicTol(this)
    micLab <- attr(mic, "label")
    div <- getDivTol(this)
    divLab <- attr(div, "label")
    psi <- getPsiTol(this)
    psiLab <- attr(psi, "label")
    
    msg2 <- paste("Convergence criteria: ",
                  "\n- ", micLab, "\t\t< ", mic, 
                  "\n- ", divLab, "\t\t< ", div, 
                  "\n- ", psiLab, "\t\t< ", psi, sep="")
    rm(psi, mic, div)
    s <- c(s, msg2, msg)
    s <- c(s, "")
  }

  ## confidence intervals
  psi <- getPsi(this)
  alpha <- 1-getConfLevel(this)
  n <- nrow(getObs(this))
  CI <- psi+c(-1, 1)*sic*qnorm(1-alpha/2)/sqrt(n)
  CI <- signif(CI, 3)
  s <- c(s, sprintf("%s-confidence interval:\t[%s, %s]", 1-alpha, CI[1], CI[2]))
  
  ## tests
  ts1 <- sqrt(n)*(psi-0)/sic
  pval1 <- 2*(1-pnorm(abs(ts1)))
  s <- c(s, sprintf("Test of \"psi(P_0)=0\":\t\tp-value = %s", signif(pval1, 3)))

  phi <- getPhi(this)
  ts2 <- sqrt(n)*(psi-phi)/getSicAlt(this)
  pval2 <- 2*(1-pnorm(abs(ts2)))
  s <- c(s, sprintf("Test of \"psi(P_0)=phi(P_0)\":\tp-value = %s",
                    signif(pval2, 3)),
         sprintf(" (estimated phi(P_0):\t%s)", signif(phi, 3)))

  class(s) <- "GenericSummary";
  s;
### A character string summarizing the content of the object. The summary contains:
### \itemize{
### \item{The sample size of the data set involved in the TMLE procedure.}
### \item{The value of the TMLE and its estimated standard error.}
### \item{A reminder of  the tuning of the stopping criteria,  and a report on
### the convergence of the TMLE procedure (see \code{\link{tmle.npvi}}).}
### \item{A confidence interval with default level of 95% (the level can be changed by using \code{\link{setConfLevel}}).}
### \item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=0}''.}
### \item{The \eqn{p}-value of the two-sided test of ``\eqn{\Psi(P_0)=\Phi(P_0)}'', with the estimated value of \eqn{\Phi(P_0)}.}
### }
}, private=TRUE)


###########################################################################/**
# @RdocMethod updateSigma2
# @alias updateSigma2
#
# @title "Updates the current estimation of parameter sigma2"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dev}{A @numeric, the "derivative" of parameter \var{sigma2}
#     at its current estimated value.}
#   \item{eps}{A @numeric, the "magnitude" of the increment to be performed.}
#   \item{\dots}{Not used.}
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{sigma2}.
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
setMethodS3("updateSigma2", "NPVI", function(this, dev, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'dev':
  dev <- Arguments$getNumeric(dev);

  eps <- getEpsilon(this);  
  sigma2 <- getSigma2(this);
  
  sigma21 <- sigma2 + eps*dev;
  ## sigma21 is positive because 'eps' is upper bounded by 1/supremum of
  ## absolute value of efficient influence curve
  setSigma2(this, sigma21)
})

###########################################################################/**
# @RdocMethod updateEfficientInfluenceCurve
# @alias updateEfficientInfluenceCurve
#
# @title "Returns the estimated efficient influence curve of parameter psi"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns a NPVI @object containing the estimated efficient influence
#  curve of parameter \var{psi}, in the form of a @matrix with 3 columns:
#  the two components of the curve, and their sum.
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
setMethodS3("updateEfficientInfluenceCurve", "NPVI", function(this, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Retrieve arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  obs <- getObs(this)
  g <- getG(this);
  mu <- getMu(this);
  theta <- getTheta(this);
  theta0 <- getTheta0(this);
  fX <- getFX(this)
  fY <- getFY(this)
  sigma2 <- getSigma2(this);
  psi <- getPsi(this);
  
  thetaXW <- theta(obs[, c("X", "W")]);
  W <- obs[, "W", drop=FALSE]
  theta0W <- theta0(W);
  muW <- mu(W);
  gW <- g(W);

  X <- fX(obs)
  Y <- fY(obs)
          
  D1 <- X * (thetaXW - theta0W - X * psi);
  D2 <- (Y - thetaXW) * (X - muW/gW*(X==0));
  verbose && summary(verbose, D1);
  verbose && summary(verbose, D2);
  verbose && print(verbose, sigma2);

  eic1 <- D1 / sigma2;
  eic2 <- D2 / sigma2;
  eic <- eic1 + eic2;
  verbose && summary(verbose, eic);

  this$.efficientInfluenceCurve <- cbind(eic1, eic2, eic);
})

###########################################################################/**
# @RdocMethod estimateEpsilon
# @alias estimateEpsilon
#
# @title "Estimates the fluctuation parameter"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{cleverCovTheta}{If @TRUE, theta is updated using the corresponding
#     "clever covariate".}
#   \item{bound}{A @numeric, upper bound for the magnitude of the change
#     in the estimated parameter}
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns the estimated flucuation parameter. 
# }
#
# @author
#
# \seealso{
#   @seemethod "updateEfficientInfluenceCurve"
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("estimateEpsilon", "NPVI", function(this, cleverCovTheta, bound=1e-1, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## Argument 'cleverCovTheta'
  cleverCovTheta <- Arguments$getLogical(cleverCovTheta);

  eic <- getEfficientInfluenceCurve(this, verbose=verbose);
  if (cleverCovTheta) {
    eic <- eic[, "eic1"]
  } else {
    eic <- eic[, "eic"];
  }
  verbose && summary(verbose, eic);
  
  if (sum(abs(eic)==Inf, na.rm=TRUE)>0) {
    throw("Infinite values in estimated efficient influence curve");
    ## eic[abs(eic)==Inf] <- NA;
  }

  theMin <- min(eic)
  theMax <- max(eic)
  if (theMin > 0) {
    interval <- c(-1/(1.001*theMax), 1e3)
  }
  if (theMax < 0) {
    interval <- c(-1e3, -1/(1.001*theMin))
  }
  if (theMin<=0 & theMax>=0) {
    interval <- c(-1/(1.001*theMax), -1/(1.001*theMin))
  }
  
  logLik <- function(epsilon) {
    sum(log(1 + epsilon * eic));
  }

  interval <- pmin(bound, pmax(-bound, interval))
  verbose && cat(verbose, "Optimization interval");
  verbose && print(verbose, interval);
  opt <- optimize(logLik, interval=interval, maximum=TRUE);

  names(opt) <- c("eps", "feps");
  eps <- opt$eps;
  attr(eps, "feps") <- opt$feps;

  eps;
})

###########################################################################/**
# @RdocMethod estimateEpsilonTheta
# @alias estimateEpsilonTheta
#
# @title "Estimates the fluctuation parameter when fluctuating 'theta'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{obs}{A @matrix of observations with 3 columns:
#     \describe{
#       \item{Y}{expression level}
#       \item{X}{DNA copy number}
#       \item{W}{DNA methylation level.}
#     }
#   }
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns the estimated flucuation parameter when fluctuating 'theta'.
# }
#
# @author
#
# \seealso{
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("estimateEpsilonTheta", "NPVI", function(this, ..., verbose=FALSE) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'obs':
  obs <- getObs(this);
  fY <- getFY(this)
  
  theta <- getTheta(this);
  H <- getHTheta(this);

  XW <- obs[, c("X", "W")]
  HXW <- H(XW)
  residuals <- fY(obs)-theta(XW)
  
  eps <- mean(residuals*HXW)/mean(HXW^2);
  feps <- eps*sum(residuals*HXW)  - (eps^2) * sum(HXW^2)/2
  attr(eps, "feps") <- feps
  
  eps
})

###########################################################################/**
# @RdocMethod updateEpsilon
# @alias updateEpsilon
#
# @title "Updates the estimate of the fluctuation parameter.
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Arguments to be passed to 'estimateEpsilon'.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns a NPVI @object containing the estimated flucuation parameter.
# }
#
# @author
#
# \seealso{
#   @seemethod "update"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("updateEpsilon", "NPVI", function(this, ..., verbose=FALSE) {
  eps <- estimateEpsilon(this, ..., verbose=verbose)

  this$.logLikIncr <- attr(eps, 'feps')
  attr(eps, 'feps') <- NULL
  this$.epsilon <- eps;
})

###########################################################################/**
# @RdocMethod updateEpsilonTheta
# @alias updateEpsilonTheta
#
# @title "Updates the estimate of the fluctuation parameter when
#   fluctuating 'theta'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{\dots}{Not used.}
#   \item{verbose}{If @TRUE, extra information is output.}
# }
#
# \value{
#  Returns a NPVI @object where the estimated flucuation parameter
#  when fluctuating 'theta' has been updated
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
setMethodS3("updateEpsilonTheta", "NPVI", function(this, ..., verbose=FALSE) {
  epsTheta <- estimateEpsilonTheta(this, ..., verbose=verbose)

  this$.logLikIncrTheta <- attr(epsTheta, 'feps')
  attr(epsTheta, 'feps') <- NULL
  this$.epsilonTheta <- epsTheta;
})


###########################################################################/**
# @RdocMethod getEpsilon
# @alias getEpsilon
#
# @title "Returns the fluctuation parameter for \acronym{TMLE}"
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
#  Returns a @numeric, the fluctuation parameter(s) for \acronym{TMLE}.
# }
#
# @author
#
# \seealso{
#   @see "update.NPVI"
#   @see "updateEpsilon.NPVI"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getEpsilon", "NPVI", function(this, ...) {
  this$.epsilon
})

###########################################################################/**
# @RdocMethod getEpsilonTheta
# @alias getEpsilonTheta
#
# @title "Returns the fluctuation parameter for \acronym{TMLE} when
#   fluctuating parameter 'theta'"
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
#  Returns a @numeric, the fluctuation parameter for \acronym{TMLE}.
# }
#
# @author
#
# \seealso{
#   @see "update.NPVI"
#   @see "updateEpsilon.NPVI"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getEpsilonTheta", "NPVI", function(this, ...) {
  this$.epsilonTheta
})


###########################################################################/**
# @RdocMethod getLogLikIncr
# @alias getLogLikIncr
#
# @title "Returns the increase in log-likelihood induced by parameter
#   updates in Targeted Maximum Likelihood Estimation"
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
#  Returns a @numeric, the increase in log-likelihood induced by parameter
#  updates in \acronym{TMLE}.
# }
#
# @author
#
# \seealso{
#   @see "update.NPVI"
#   @see "updateEfficientInfluenceCurve.NPVI"
#   @see "updateEpsilon.NPVI"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getLogLikIncr", "NPVI", function(this, ...) {
  this$.logLikIncr
})


###########################################################################/**
# @RdocMethod getLogLikIncrTheta
# @alias getLogLikIncrTheta
#
# @title "Returns the increase in log-likelihood induced by parameter
#   updates in Targeted Maximum Likelihood Estimation"
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
#  Returns a @numeric, the increase in log-likelihood induced when
#  fluctuating \var{theta} in \acronym{TMLE}.
# }
#
# @author
#
# \seealso{
#   @see "update.NPVI"
#   @see "updateEfficientInfluenceCurve.NPVI"
#   @see "updateEpsilon.NPVI"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getLogLikIncrTheta", "NPVI", function(this, ...) {
  this$.logLikIncrTheta
})


###########################################################################/**
# @RdocMethod getEfficientInfluenceCurve
# @alias getEfficientInfluenceCurve
#
# @title "Returns the estimated value of the efficient influence curve
#   that has been used for parameter updates by Targeted Maximum Likelihood
#   Estimation"
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
#  Returns a @numeric, the estimated value of the efficient influence curve
#   that has been used for parameter updates.
# }
#
# @author
#
# \seealso{
#   @see "update.NPVI"
#   @see "updateEfficientInfluenceCurve.NPVI"
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getEfficientInfluenceCurve", "NPVI", function(this, ...) {
  this$.efficientInfluenceCurve;
})

###########################################################################/**
# @RdocMethod getWeightsW
# @alias getWeightsW
#
# @title "Returns the value of parameter 'weightsW'."
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
#  Returns a @vector, the value of parameter \var{weightsW}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
#*/###########################################################################
setMethodS3("getWeightsW", "NPVI", function(this, ...) {
  this$.weightsW
})



###########################################################################/**
# @RdocMethod setWeightsW
# @alias setWeightsW
#
# @title "Sets the value of parameter 'weightsW'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mu}{A @vector of non-negative weights (it is not necessary that it sums to 1). }
# }
#
# \value{  Returns  a NPVI  @object  containing  the  updated vector of
#         weights for marginal simulation of W. }
#
# @author
#
# \seealso{
#   @seeclass
#   @seeMethod getWeightsW
# }
#
#*/###########################################################################
setMethodS3("setWeightsW", "NPVI", function(this, weightsW, ...) {
  ## Argument 'weightsW':
  weightsW <- Arguments$getNumerics(weightsW, range=c(0, Inf))
  nr <- nrow(getObs(this))
  if (length(weightsW) != nr) {
    throw("Length of 'weightsW' must match number of observations!")
  }
  this$.weightsW <- weightsW
})

###########################################################################/**
# @RdocMethod updateWeightsW
# @alias updateWeightsW
#
# @title "Updates the current value of parameter 'weightsW'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{effICW}{A @function, the conditional expectation given W of
#      projection of efficient influence curve on functions of (W,X). }
# }
#
# \value{
#  Returns a NPVI @object containing the updated \var{weightsW}.
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
setMethodS3("updateWeightsW", "NPVI", function(this, effICW, ...) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'effICW':
  if (!is.null(effICW) & (mode(effICW)!="function")) {
    throw("Argument \var{effICW} should be of mode 'function', not ", mode(effICW));
  }

  if (!is.null(effICW)) {
    obs <- getObs(this)
    eps <- getEpsilon(this);  
    weightsW <- getWeightsW(this)*(1+eps*effICW(obs[, "W"]))

    setWeightsW(this, weightsW)
  }
})

setMethodS3("getConv", "NPVI", function(this, ...) {
  this$.conv;
})

setMethodS3("updateConv", "NPVI", function(x, B, ...) {
  this <- x;  ## To please R CMD check

  n <- nrow(getObs(this))
  mic.tol <- getMicTol(this)
  div.tol <- getDivTol(this)
  psi.tol <- getPsiTol(this)
  ## extracting the relevant components of history
  ## (whether 'cleverCovTheta' is TRUE or FALSE)
  hist <- getHistory(this)
  step <- as.integer(gsub("step([0-9]+)", "\\1", rownames(hist)))
  idxs <- c(which(diff(step)==1), length(step))
  hist <- tail(hist[idxs, ], 2)

  conv <- FALSE
  msg <- NULL
  
  ## Assessing values of scaled 'mic'
  if (!is.na(mic.tol)) {
    mic <- hist[2, "mic"]
    Zmic <- hist[2, "psi.sd"]/sqrt(n)
    scaledMic <- mic/Zmic
    micLab <- attr(mic.tol, "label")
    if (abs(scaledMic) < mic.tol) {
      conv <- TRUE;
      msg1 <- paste(micLab, " is within ", mic.tol, "-tolerance", sep="");
      msg <- c(msg, msg1)
    }
  }
  
  ## Testing value of 'div'
  if (!is.na(div.tol)) {
    div <- hist[2, "div"]
    divLab <- attr(div.tol, "label")
    if (!is.na(div)) {
      if (div < div.tol) {
        conv <- TRUE;
        msg1 <- paste(divLab, " is within ", div.tol, "-tolerance", sep="");
        msg <- c(msg, msg1)
      }
    }
  }

  ## Comparing successive values of 'psi'
  if (!is.na(psi.tol)) {
    psiLab <- attr(psi.tol, "label")
    psi <- hist[, "psi"]
    psi.sd <- hist[, "psi.sd"]
    crit <- abs(diff(psi))/sqrt(sum(psi.sd^2))*sqrt(B)
    if (crit < psi.tol) {
      conv <- TRUE;
      msg1 <- paste(psiLab, " are within ", psi.tol, "-tolerance", sep="");
      msg <- c(msg, msg1)
    }
  }

  if (!conv) {
    msg <- "Convergence not reached (yet)"
  }

  
  attr(conv, "msg") <- msg
  this$.conv <- conv
})



############################################################################
## HISTORY:
## 2011-05-18
## o Added standard deviation of EIC and divergence to history
## o Added functions 'setDivergence' and 'getDivergence'
## 2011-05-02
## o Made distinct files for functions related to g, mu and theta.
## o Added parameter 'weightsW'
## 2011-04-22
## o Moved function 'updatePsi' to its own file.
## o Introduced 'tabulate' parameter...
## o Call to new function 'validateArgumentObs'
## 2011-03-21
## o Added 'family' parameter and related functions 'getFamily', 'setFamily',
##   'getObs', and modified accordingly functions 'getG' and 'setG' etc.
## o Added methods 'getGW', 'getMuW', and 'getThetaXW'.
## o Removed 'obs' and 'family' from arguments of 'updatePsi'
## 2011-02-28
## o Added element 'step' and corresponding method to structure.
## 2011-02-17
## o Added element 'history' and corresponding methods.
## o Added methods 'updatePsi', 'setG', 'setMu', 'setTheta'.
## o Major updates to class definition.
## 2011-02-16
## o Merged classes NPVI and TMLE.NPVI.
## o Removed 'paramList' from class NPVI.
## o Renamed EstimateECM into NPVI.
## 2011-02-08
## o Added parameters 'gmax', 'mumin', 'mumax', 'thetamin', 'thetamax',
## the related functions, and adapted code accordingly.
## 2011-02-04
## o Added 'getPsiPn'
## 2011-01-31
## o Added 'logLikIncr'
## 2011-01-24
## o Removed parameter 'theta0'
## 2011-01-20
## o modified method 'estimateEpsilon' (taking argument 'cleverCovTheta' into account)
## o added method 'estimateEpsilonTheta'
## o added function and method related to 'HTheta'
## 2010-12-31
## o added thresholding procedure of 'g'
## 2010-11-29
## o 'estimateInfluenceCurve' now returns a matrix with the components
##   D1, D2 and D1+D2.
## 2010-08-03
## o Added parameter 'theta0'.
## 2010-07-07
## o Created.
############################################################################

