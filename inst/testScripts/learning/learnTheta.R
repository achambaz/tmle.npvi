learnTheta <- function(obs, light=TRUE, ...) {
  varNames <- setdiff(colnames(obs), "Y")
  theFormula <- paste(varNames, collapse="*")
  theFormula <- paste("Y ~", theFormula, sep=" ")
  ## formula <- as.formula(Y~X*W);
  formula <- as.formula(theFormula);
  
  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(XW) {
    predict(fit, newdata=as.data.frame(XW), type="response");
  }
  attr(foo, 'fit') <- fit;

  foo
}
