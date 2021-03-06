SL.nnls <- function#h2oEnsembleLearning Wrapper For Nonnegative Least Squares Metalearning
### SuperLearner-API wrapper  function to  support NNLS for  metalearning when
### using h2oEnsembleLearning
(Y,
### The outcome in the training data set. Must be a numeric \code{vector}.
 X,
### The   predictor  variables   in   the  training   data   set,  usually   a
### \code{data.frame}.
 newX,
### The predictor variables  in the validation data set.  The structure should
### match \code{X}. If missing, uses \code{X} for \code{newX}.
 family,
### Currently  allows   "gaussian"  or   "binomial"  to  describe   the  error
### distribution.
 obsWeights,
###  Optional observation weights variable. It is passed to the prediction and
###  screening algorithms, but many of  the built in wrappers ignore (or can't
###  use) the  information.  If you  are using observation weights,  make sure
###  the library you specify uses the information.
 ...
### Not used.
 ) {
  ##seealso<< predict.SL.nnls, method.NNLS
  requireNamespace("nnls")
  fit.nnls <- nnls::nnls(sqrt(obsWeights)*as.matrix(X), sqrt(obsWeights)*Y) 
  initCoef <- coef(fit.nnls)
  initCoef[is.na(initCoef)] <- 0
  if (sum(initCoef) > 0) {
    coef <- initCoef/sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call.=FALSE)
    coef <- initCoef
  }
  pred <- crossprod(t(as.matrix(newX)), coef)
  fit <- list(object=fit.nnls)
  class(fit) <- "SL.nnls"
  out <- list(pred=pred, fit=fit)
  return(out)
### A fitted object.
}

predict.SL.nnls <- function#h2oEnsembleLearning Wrapper For Nonnegative Least Squares Metalearning
### h2oEnsembleLearning wrapper for nonnegative least squares metalearning
(object,
### A fitted object as given by \code{SL.nnls}.
 newdata,
### The predictor variables for which predictions are wished.
 ...
### Not used.
 ) {
  ##seealso<< SL.nnls, method.NNLS
  initCoef <- coef(object$object)
  initCoef[is.na(initCoef)] <- 0
  if (sum(initCoef) > 0) {
    coef <- initCoef/sum(initCoef)
  } else {
    warning("All algorithms have zero weight", call.=FALSE)
    coef <- initCoef
  }
  pred <- crossprod(t(as.matrix(newdata)), coef)
  return(pred)
### Returns a vector of predictions for \code{newdata} derived from the fitted
### object \code{object}.
}


h2o.randomForest.1000x100 <- function#Builds A Random Forest Model On An H2OFrame
### Builds a random forest model  on an H2OFrame, with parameter \code{ntrees}
### and \code{nbins} set to \code{1000} and \code{100}.
(...,
### All parameters of \code{h2o.randomForest} except \code{ntrees} and \code{nbins}.
 ntrees=1000,
### Parameter \code{ntrees} of \code{h2o.randomForest}, set to \code{1000}.
 nbins=100
### Parameter \code{nbins} of \code{h2o.randomForest}, set to \code{100}.
 ){
  ##seealso<< h2o.randomForest
  h2oEnsemble::h2o.randomForest.wrapper(..., ntrees=ntrees, nbins=nbins, seed=1)
### Creates a H2OModel object of the right type.  
}

h2o.glm.alpha.00 <- function#H2O Generalized Linear Models
### Fit  a  generalized  linear  model,  with parameter  \code{alpha}  set  to
### \code{0.0}.
(...,
### All parameters of \code{h2o.glm} except \code{alpha}.
 alpha=0.0
### Parameter \code{alpha} of \code{h2o.glm}, set to \code{0.0}.
 ){
  ##seealso<< h2o.glm
  h2oEnsemble::h2o.glm.wrapper(..., alpha=alpha)
### Creates a H2OModel object of the right type.  
}

h2o.glm.alpha.05 <- function#H2O Generalized Linear Models
### Fit  a  generalized  linear  model,  with parameter  \code{alpha}  set  to
### \code{0.5}.
(...,
### All parameters of \code{h2o.glm} except \code{alpha}.
 alpha=0.5
### Parameter \code{alpha} of \code{h2o.glm}, set to \code{0.5}.
 ){
  ##seealso<< h2o.glm
  h2oEnsemble::h2o.glm.wrapper(..., alpha=alpha)
### Creates a H2OModel object of the right type.  
}

h2o.glm.alpha.10 <- function#H2O Generalized Linear Models
### Fit  a  generalized  linear  model,  with parameter  \code{alpha}  set  to
### \code{1.0}.
(...,
### All parameters of \code{h2o.glm} except \code{alpha}.
 alpha=1.0
### Parameter \code{alpha} of \code{h2o.glm}, set to \code{1.0}.
 ){
  ##seealso<< h2o.glm
  h2oEnsemble::h2o.glm.wrapper(..., alpha=alpha)
### Creates a H2OModel object of the right type.  
}

h2o.deeplearning.Rectifier <- function#Deep Learning Neural Network
### Performs  Deep Learning neural  networks on  an H2OFrame,  with parameters
### \code{activation}] and  \code{hidden} set to  "Rectifier" and \code{c(500,
### 500)}.
(...,
### All  parameters   of  \code{h2o.deeplearning}  except   \code{hidden}  and
### \code{activation}.
 hidden=c(500, 500),
### Parameter  \code{hidden} of  \code{h2o.deeplearning}, set  to \code{c(500,
### 500)}.
 activation="Rectifier"
### Parameter \code{activation} of \code{h2o.deeplearning}, set to "Rectifier".
 ) {
  ##seealso<< h2o.deeplearning
  h2oEnsemble::h2o.deeplearning.wrapper(..., hidden=hidden, activation=activation,
                                        seed=1)
### Creates a H2OModel object of the right type.  
}

h2o.deeplearning.Tanh <- function#Deep Learning Neural Network
### Performs  Deep Learning neural  networks on  an H2OFrame,  with parameters
### \code{activation}] and  \code{hidden} set to "Tanh"  and \code{c(200, 200,
### 200)}.
(...,
### All  parameters   of  \code{h2o.deeplearning}  except   \code{hidden}  and
### \code{activation}.
 hidden=c(200, 200, 200),
### Parameter  \code{hidden} of  \code{h2o.deeplearning}, set  to \code{c(200,
### 200, 200)}.
 activation="Tanh"
### Parameter \code{activation} of \code{h2o.deeplearning}, set to "Rectifier".
 ) {
  ##seealso<< h2o.deeplearning
  h2oEnsemble::h2o.deeplearning.wrapper(..., hidden=hidden, activation=activation,
                                        seed=1)
### Creates a H2OModel object of the right type.  
}
