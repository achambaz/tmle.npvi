SL.DSA <- function
### Prediction algorithm wrapper for SuperLearner.
(Y, X, newX, family, obsWeights, maxsize = ncol(X),
 maxorderint = 1, maxsumofpow = 1, Dmove = TRUE, Smove = TRUE,
 vfold = 5, ...) {
  ##seealso<< SL.DSA.2, predict.SL.DSA
  tryCatch(require(DSA), warning = function(...) {
    stop("you have selected DSA as a library algorithm but do not have the DSA package installed")
  })
  dsaweights <- matrix(obsWeights, nrow = (vfold + 1), ncol = nrow(X),
                       byrow = TRUE)
  fit.DSA <- DSA(Y ~ 1, data = data.frame(Y, X),
                 family = family, maxsize = maxsize, maxorderint = maxorderint,
                 maxsumofpow = maxsumofpow, Dmove = Dmove, Smove = Smove,
                 vfold = vfold, weights = dsaweights)
  pred <- predict(fit.DSA, newdata = newX)
  if (family$family == "binomial") {
    pred <- 1/(1 + exp(-pred))
  }
  fit <- list(object = fit.DSA)
  class(fit) <- "SL.DSA"
  foo <- list(pred = pred, fit = fit)
  return(foo)
### Returns a fitted object.
}
environment(SL.DSA) <- asNamespace("SuperLearner")


SL.DSA.2 <- function
### Prediction algorithm wrapper for SuperLearner.
(..., X, maxsize = 2 * ncol(X), maxorderint = 2,
 maxsumofpow = 2, Smove = FALSE, vfold = 10) {
  ##seealso<< SL.DSA, predict.SL.DSA
  out <- SL.DSA(..., X = X, maxsize = maxsize, maxorderint = maxorderint,
                maxsumofpow = maxsumofpow, Smove = Smove, vfold = vfold)
  return(out)
### Returns a fitted object.
}
environment(SL.DSA.2) <- asNamespace("SuperLearner")



predict.SL.DSA <- function
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, family, X = NULL, Y = NULL, ...) {
  ##seealso<< SL.DSA, SL.DSA.2
  tryCatch(require(DSA), warning = function(...) {
    stop("you have selected DSA as a library algorithm but do not have the DSA package installed")
  })
  pred <- predict(object = object$object, newdata = newdata)
  if (family$family == "binomial") {
    pred <- 1/(1 + exp(-pred))
  }
  return(pred)
### Returns a prediction.
}
environment(predict.SL.DSA) <- asNamespace("SuperLearner")


SL.glm.condExpX2givenW <- function
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnCondExpX2givenW, predict.SL.glm.condExpX2givenW
  varNames <- names(X)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep=" ")
  ## formula.glm.condExpX2givenW <- as.formula(Y~W+I(W^2));
  formula.glm.condExpX2givenW <- as.formula(theFormula);
  
  fit.glm <- glm(formula.glm.condExpX2givenW, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.condExpX2givenW")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
environment(SL.glm.condExpX2givenW) <- asNamespace("SuperLearner")


predict.SL.glm.condExpX2givenW <- function
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
(object, newdata, ...) {
  ##seealso<< SL.glm.condExpX2givenW
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}
environment(predict.SL.glm.condExpX2givenW) <- asNamespace("SuperLearner")


SL.glm.condExpXYgivenW <- function
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{XY} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnCondExpXYgivenW, predict.SL.glm.condExpXYgivenW
  varNames <- names(X)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep=" ")
  ## formula.glm.condExpXYgivenW <- as.formula(Y~W+I(W^2));
  formula.glm.condExpXYgivenW <- as.formula(theFormula);
  
  fit.glm <- glm(formula.glm.condExpXYgivenW, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.condExpXYgivenW")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
environment(SL.glm.condExpXYgivenW) <- asNamespace("SuperLearner")


predict.SL.glm.condExpXYgivenW <- function
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.condExpXYgivenW
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}
environment(predict.SL.glm.condExpXYgivenW) <- asNamespace("SuperLearner")


SL.glm.g <- function
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional probability of \eqn{X=0} given \eqn{W}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnG, predict.SL.glm.g
  varNames <- names(X)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep=" ")
  ## formula.glm.g <- as.formula(Y~W+I(W^2));
  formula.glm.g <- as.formula(theFormula);
  fit.glm <- glm(formula.glm.g, data = X, family = family,
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.g")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
environment(SL.glm.g) <- asNamespace("SuperLearner")


predict.SL.glm.g <- function
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.g
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    return(out)
### Returns a prediction.
}
environment(predict.SL.glm.g) <- asNamespace("SuperLearner")


SL.glm.theta <- function
### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{Y} given \eqn{(X,W)}.
(Y, X, newX, family, obsWeights, ...) {
  ##seealso<< learnTheta, predict.SL.glm.theta
  varNames <- names(X)
  theFormula <- paste(varNames, collapse="*")
  theFormula <- paste("Y ~", theFormula, sep=" ")
  ## formula.glm.theta <- as.formula(Y~X*W);
  formula.glm.theta <- as.formula(theFormula);
  
  fit.glm <- glm(formula.glm.theta, data = X, family = family, 
                 weights = obsWeights)
  pred <- predict(fit.glm, newdata = newX, type = "response")
  fit <- list(object = fit.glm)
  class(fit) <- c("SL.glm.theta")
  out <- list(pred = pred, fit = fit)
  return(out)
### Returns a fitted object.
}
environment(SL.glm.theta) <- asNamespace("SuperLearner")



predict.SL.glm.theta <- function
### Prediction algorithm wrapper for SuperLearner.
(object, newdata, ...) {
  ##seealso<< SL.glm.theta
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    return(out)
### Returns a prediction.
}
environment(predict.SL.glm.theta) <- asNamespace("SuperLearner")


## -----------------------------------------------------------------------
## setting the different libraries to use when 'flavor' is "superLearning"
## -----------------------------------------------------------------------

## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA")
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");

SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm")

library(SuperLearner)
library(randomForest)
library(e1071)
## library(polspline)
## library(DSA)
## library(glmnet)

learnTheta.library <- c("SL.glm.theta", "SL.polymars", SL.library);

learnG.library <- c("SL.glm.g", SL.library);

learnMuAux.library <- c(SL.library);

learnDevG.library <- c(SL.library)

learnDevMu.library <- c(SL.library)

learnDevTheta.library <- c(SL.library)

learnCondExpXYgivenW.library <- c("SL.glm.condExpXYgivenW", SL.library)

learnCondExpX2givenW.library <- c("SL.glm.condExpX2givenW", SL.library)


### List of default libraries of algorithms to use in \code{tmle.npvi} when \code{flavor} is set to "learning".
superLearningLib <- list(learnCondExpX2givenW=learnCondExpX2givenW.library,
                         learnCondExpXYgivenW=learnCondExpXYgivenW.library,
                         learnDevG=learnDevG.library,
                         learnDevMu=learnDevMu.library,
                         learnDevTheta=learnDevTheta.library,
                         learnG=learnG.library,
                         learnMuAux=learnMuAux.library,
                         learnTheta=learnTheta.library)

### Default library of algorithms to use in \code{tmle.npvi} when argument \code{flavor} is set to "superLearning". 
SL.library <- unique(c(learnTheta.library, learnG.library, learnMuAux.library,
                       learnDevTheta.library, learnDevG.library, learnDevMu.library,
                       learnCondExpXYgivenW.library, learnCondExpX2givenW.library))

rm(learnTheta.library,           learnG.library,           learnMuAux.library,
   learnDevTheta.library,          learnDevG.library,         learnDevMu.library,
   learnCondExpXYgivenW.library, learnCondExpX2givenW.library)


