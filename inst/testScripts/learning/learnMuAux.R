learnMuAux <- function(obs, light=TRUE, ...) {
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("X ~", theFormula, "+", theFormula2, sep=" ")
  ## formula <- as.formula(X~W+I(W^2));
  formula <- as.formula(theFormula)

  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;
  foo
}
