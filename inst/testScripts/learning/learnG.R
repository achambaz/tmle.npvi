learnG <- function(obs, theX0=0, light=TRUE, ...) {
  ##
  ## version glm
  ##
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  theFormula <- paste(varNames, collapse=" + ")
  theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
  theFormula <- paste("I(X==theX0) ~", theFormula, "+", theFormula2, sep=" ")
  ## formula <- as.formula(I(X==theX0)~W);
  formula <- as.formula(theFormula)
    
  fit <- glm(formula, data=as.data.frame(obs), family="binomial");
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    ## predict(fit, newdata=data.frame(W=W), type="response")
    predict(fit, newdata=data.frame(W), type="response")
  }
  attr(foo, 'fit') <- fit;

  foo
}

## environment(learnG) <- asNamespace("TMLE.NPVI")
