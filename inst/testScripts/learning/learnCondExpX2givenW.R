learnCondExpX2givenW <- function(obs, light=TRUE) {
  ## W <- obs[, "W"];
  W <- extractW(obs)
  varNames <- colnames(W)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("I(X^2) ~", theFormula, "+", theFormula2, sep="")
  ## formula <- as.formula(I(X^2)~W+I(W^2));
  formula <- as.formula(theFormula)
  rm(W);

  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;
  foo;
}

