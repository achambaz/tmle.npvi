learnDevG <- function(obs, effIC1, gW, light=TRUE, verbose=FALSE, ...) {
  ## W <- obs[, "W"];
  W <- extractW(obs)
  X <- obs[, "X"];
  Z <- effIC1 * ( (X==0) - gW );

  obsZ <- cbind(obs, Z=Z);
  verbose && str(verbose, obsZ);

  varNames <- colnames(W)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("Z ~", theFormula, "+", theFormula2, sep="")
  ## formula <- as.formula(Z~W+I(W^2));
  formula <- as.formula(theFormula)

  fit <- glm(formula, data=as.data.frame(obsZ), family=gaussian);
  rm(W, X, Z, obsZ);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;

  foo
}
