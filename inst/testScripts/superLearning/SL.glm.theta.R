SL.glm.theta <- function(Y, X, newX, family, obsWeights, ...) {
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
}
environment(SL.glm.theta) <- asNamespace("SuperLearner")


predict.SL.glm.theta <- function (object, newdata, ...) {
    out <- predict(object = object$object, newdata = newdata, 
        type = "response")
    out
}
environment(predict.SL.glm.theta) <- asNamespace("SuperLearner")
