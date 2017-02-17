#' @importFrom stats as.formula glm
#' @importFrom SuperLearner SL.glm SL.randomForest SL.gam SL.polymars
SL.glm.condExpX2givenW <- function(Y, X, newX, family, obsWeights, ...) {
    varNames <- names(X)
    if (length(varNames)>20) {
        varNames <- varNames[1:20]
        warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.condExpX2givenW'"))
    }
    theFormula <- paste(varNames, collapse=" + ")
    if (length(varNames)<=10) {
        theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
        theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep="")
    } else {
        theFormula <- paste("Y ~", theFormula, sep="")
    }
    ## formula.glm.condExpX2givenW <- as.formula(Y~W+I(W^2));
    formula.glm.condExpX2givenW <- as.formula(theFormula);

    fit.glm <- glm(formula.glm.condExpX2givenW, data = X, family = family,
                   weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- c("SL.glm.condExpX2givenW")
    out <- list(pred = pred, fit = fit)
    return(out)
}



predict.SL.glm.condExpX2givenW <- function(object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata,
                   type = "response")
    return(out)
}


SL.glm.condExpXYgivenW <- function(Y, X, newX, family, obsWeights, ...) {
    varNames <- names(X)
    if (length(varNames)>20) {
        varNames <- varNames[1:20]
        warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.condExpXYgivenW'", collapse=""))
    }
    theFormula <- paste(varNames, collapse=" + ")
    if (length(varNames)<=10) {
        theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
        theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep="")
    } else {
        theFormula <- paste("Y ~", theFormula, sep="")
    }
    ## formula.glm.condExpXYgivenW <- as.formula(Y~W+I(W^2));
    formula.glm.condExpXYgivenW <- as.formula(theFormula);

    fit.glm <- glm(formula.glm.condExpXYgivenW, data = X, family = family,
                   weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- c("SL.glm.condExpXYgivenW")
    out <- list(pred = pred, fit = fit)
    return(out)
}


predict.SL.glm.condExpXYgivenW <- function(object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata,
                   type = "response")
    return(out)
}



SL.glm.g <- function(Y, X, newX, family, obsWeights, ...) {
    varNames <- names(X)
    if (length(varNames)>20) {
        varNames <- varNames[1:20]
        warning(paste("Using only", paste(varNames, collapse=", "), "in 'SL.glm.g'"))
    }
    theFormula <- paste(varNames, collapse=" + ")
    if (length(varNames)<=10) {
        theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
        theFormula <- paste("Y ~", theFormula, "+", theFormula2, sep=" ")
    } else {
        theFormula <- paste("Y ~", theFormula, sep="")
    }
    ## formula.glm.g <- as.formula(Y~W+I(W^2));
    formula.glm.g <- as.formula(theFormula);



    fit.glm <- glm(formula.glm.g, data = X, family = family,
                   weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- c("SL.glm.g")
    out <- list(pred = pred, fit = fit)
    return(out)
}

#' @export
predict.SL.glm.g <- function(object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata,
                   type = "response")
    return(out)
}


SL.glm.theta <- function(Y, X, newX, family, obsWeights, ...) {
    varNames <- names(X)
    if (length(varNames)>20) {
        varNames <- varNames[1:20]
        warning(paste("Using only 'X' and", paste(varNames, collapse=", "), "in 'SL.glm.theta'"))
    }
    theFormula <- paste(varNames, collapse=" + ")
    if (length(varNames)<=10) {
        theFormula <- paste("Y ~ X+", theFormula, "+ X*(",
                            theFormula, ")", sep="")
    } else {
        theFormula <- paste("Y ~ X+", theFormula, sep="")
    }
    ## formula.glm.theta <- as.formula(Y~X*W);
    formula.glm.theta <- as.formula(theFormula);

    fit.glm <- glm(formula.glm.theta, data = X, family = family,
                   weights = obsWeights)
    pred <- predict(fit.glm, newdata = newX, type = "response")
    fit <- list(object = fit.glm)
    class(fit) <- c("SL.glm.theta")
    out <- list(pred = pred, fit = fit)
    return(out)
}
## fyi:
## environment(SL.glm.theta) <- asNamespace("SuperLearner")
## environment(predict.SL.glm.theta) <- asNamespace("SuperLearner")


predict.SL.glm.theta <- function(object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata,
                   type = "response")
    return(out)
}
