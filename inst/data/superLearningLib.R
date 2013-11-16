### Prediction algorithm wrapper for SuperLearner
SL.DSA <- function
(Y, X, newX, family, obsWeights, maxsize = ncol(X),
 maxorderint = 1, maxsumofpow = 1, Dmove = TRUE, Smove = TRUE,
 vfold = 5, ...) {
  ###seealso<< SL.DSA.2, predict.SL.DSA
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

### Prediction algorithm wrapper for SuperLearner
SL.DSA.2 <- function
(..., X, maxsize = 2 * ncol(X), maxorderint = 2,
 maxsumofpow = 2, Smove = FALSE, vfold = 10) {
###seealso<< SL.DSA, predict.SL.DSA
  out <- SL.DSA(..., X = X, maxsize = maxsize, maxorderint = maxorderint,
                maxsumofpow = maxsumofpow, Smove = Smove, vfold = vfold)
  return(out)
### Returns a fitted object.
}
environment(SL.DSA.2) <- asNamespace("SuperLearner")


### Prediction algorithm wrapper for SuperLearner
predict.SL.DSA <- function
(object, newdata, family, X = NULL, Y = NULL, ...) {
  tryCatch(require(DSA), warning = function(...) {
###seealso<< SL.DSA, SL.DSA.2
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


### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
SL.glm.condExpX2givenW <- function
(Y, X, newX, family, obsWeights, ...) {
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


### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{X^2} given \eqn{W}.
predict.SL.glm.condExpX2givenW <- function
(object, newdata, ...) {
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}
environment(predict.SL.glm.condExpX2givenW) <- asNamespace("SuperLearner")


### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{XY} given \eqn{W}.
SL.glm.condExpXYgivenW <- function
(Y, X, newX, family, obsWeights, ...) {
###seealso predict.SL.glm.condExpXYgivenW
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

### Prediction algorithm wrapper for SuperLearner
predict.SL.glm.condExpXYgivenW <- function
(object, newdata, ...) {
###seealso SL.glm.condExpXYgivenW
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  return(out)
### Returns a prediction.
}
environment(predict.SL.glm.condExpXYgivenW) <- asNamespace("SuperLearner")


### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional probability of \eqn{X=0} given \eqn{W}.
SL.glm.g <- function
(Y, X, newX, family, obsWeights, ...) {
###seealso predict.SL.glm.g
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

### Prediction algorithm wrapper for SuperLearner
predict.SL.glm.g <- function
(object, newdata, ...) {
  ###seealso SL.glm.g
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    return(out)
### Returns a prediction.
}
environment(predict.SL.glm.g) <- asNamespace("SuperLearner")


### Prediction algorithm  wrapper for SuperLearner, for the  estimation of the
### conditional expectation of \eqn{Y} given \eqn{(X,W)}.
SL.glm.theta <- function
(Y, X, newX, family, obsWeights, ...) {
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


### Prediction algorithm wrapper for SuperLearner
predict.SL.glm.theta <- function
(object, newdata, ...) {
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

learnTheta <- c("SL.glm.theta", "SL.polymars", SL.library);

learnG <- c("SL.glm.g", SL.library);

learnMuAux <- c(SL.library);

learnDevG <- c(SL.library)

learnDevMu <- c(SL.library)

learnDevTheta <- c(SL.library)

learnCondExpXYgivenW <- c("SL.glm.condExpXYgivenW", SL.library)

learnCondExpX2givenW <- c("SL.glm.condExpX2givenW", SL.library)

SL.library <- unique(c(learnTheta, learnG, learnMuAux,
                       learnDevTheta, learnDevG, learnDevMu,
                       learnCondExpXYgivenW, learnCondExpX2givenW))
