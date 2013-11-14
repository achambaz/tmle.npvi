SL.glm.g <- function(Y, X, newX, family, obsWeights, ...) {
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
}
environment(SL.glm.g) <- asNamespace("SuperLearner")

predict.SL.glm.g <- function (object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    out
}
environment(predict.SL.glm.g) <- asNamespace("SuperLearner")
