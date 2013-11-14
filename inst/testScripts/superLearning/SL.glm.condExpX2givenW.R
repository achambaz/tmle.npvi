SL.glm.condExpX2givenW <- function(Y, X, newX, family, obsWeights, ...) {
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
}
environment(SL.glm.condExpX2givenW) <- asNamespace("SuperLearner")


predict.SL.glm.condExpX2givenW <- function (object, newdata, ...) {
  out <- predict(object = object$object, newdata = newdata, 
                 type = "response")
  out
}
environment(predict.SL.glm.condExpX2givenW) <- asNamespace("SuperLearner")
