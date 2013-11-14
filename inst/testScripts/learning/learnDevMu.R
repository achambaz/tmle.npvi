learnDevMu <- function(obs, effIC1, muW, light=TRUE, verbose=FALSE, ...) {
  ## W <- obs[, "W"];
  W <- extractW(obs)
  Z <- (obs[, "X"]-muW) * effIC1;

  obsZ <- cbind(obs, Z=Z);
  verbose && str(verbose, obsZ);

  varNames <- colnames(W)
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("Z ~", theFormula, "+", theFormula2, sep="")
  ## formula <- as.formula(Z~W+I(W^2));
  formula <- as.formula(theFormula)

  fit <- glm(formula, data=as.data.frame(obsZ), family=gaussian);
  rm(W, Z, obsZ);
  if (light) {
    fit <- getLightFit(fit);
  }
  
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;

  foo
}
