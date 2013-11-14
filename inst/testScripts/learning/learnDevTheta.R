learnDevTheta <- function(obs, thetaXW, light=TRUE, verbose=FALSE, ...) {
  obsZ <- cbind(obs, Z=thetaXW);
  verbose && str(verbose, obsZ);

  varNames <- setdiff(colnames(obs), "Y")
  theFormula <- paste(varNames, collapse="*")
  theFormula <- paste("I((Y-Z)^2) ~", theFormula, sep=" ")
  ## formula <- as.formula(I((Y-Z)^2)~X*W);
  formula <- as.formula(theFormula)
  
  ## family <- Gamma(link="log");
  family <- gaussian();
  fit <- glm(formula, data=as.data.frame(obsZ), family=family);
  rm(obsZ)
  if (light) {
    fit <- getLightFit(fit);
  }
  
  foo <- function(XW) {
    predict(fit, newdata=as.data.frame(XW), type="response");
  }

  attr(foo, 'fit') <- fit;
  foo
}
