learnCondExpX2givenW <- function#Estimation of Cond. Expect. of X^2 Given W
### Function for  the estimation of  the conditional expectation  of \eqn{X^2}
### given \eqn{W} when \code{flavor} is set to "learning".
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 light=TRUE
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
 ) {
  ##seealso<< learnG, learnMuAux, learnTheta, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnCondExpX2givenW'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("I(X^2) ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("I(X^2) ~", theFormula, sep="")
  } 
  formula <- as.formula(theFormula)
  ## formula <- as.formula(I(X^2)~W+I(W^2));
  
  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;
  return(foo)
### Returns the fitted object.
}


learnCondExpXYgivenW <- function#Estimation of Cond. Expect. of XY Given W
### Function for  the estimation of  the conditional expectation  of \eqn{XY}
### given \eqn{W} when \code{flavor} is set to "learning".
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 light=TRUE
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
) {
  ##seealso<< learnG, learnMuAux, learnTheta, learnCondExpX2givenW, learnDevG, learnDevMu, learnDevTheta
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnCondExpXYgivenW'", collapse=""))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("I(X*Y) ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("I(X*Y) ~", theFormula, sep="")
  } 
  formula <- as.formula(theFormula)  
  ## formula <- as.formula(I(X*Y)~W+I(W^2));

  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;
  return(foo)
### Returns the fitted object.
}

learnDevG <- function#Estimation of Cond. Expect. of ((X==0)-gW)*effIC1 Given W 
### Function   for  the   estimation   of  the   conditional  expectation   of
### \code{((X==0)-gW)*effIC1}  given  \eqn{W}  when  \code{flavor} is  set  to
### "learning".
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 effIC1,
### The \code{vector}  \code{effIC1} of the  first component of  the efficient
### influence curve, as currently estimated, evaluated at our observations.
 gW,
### The \code{vector} \code{gW} of  the conditional probability that \eqn{X=0}
### given \eqn{W}, as currently estimated, evaluated at our observations.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
verbose=FALSE,
### Prescribes the amount of information  output by the function.  Defaults to
### \code{FALSE}.
 ...
### Additional arguments possibly needed.
 ) {
  ##seealso<< learnG, learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevMu, learnDevTheta
  X <- obs[, "X"];
  Z <- effIC1 * ( (X==0) - gW );

  obsZ <- cbind(obs, Z=Z);
  verbose && str(verbose, obsZ);

  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnDevG'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("Z ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("Z ~", theFormula, sep="")
  } 
  formula <- as.formula(theFormula)
  ## formula <- as.formula(Z~W+I(W^2));

  fit <- glm(formula, data=as.data.frame(obsZ), family=gaussian);
  rm(X, Z, obsZ);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;

  return(foo)
### Returns the fitted object.
}

learnDevMu <- function#Estimation of Cond. Expect. of (X-muW)*effIC1 Given W 
### Function   for  the   estimation   of  the   conditional  expectation   of
### \code{(X-muW)*effIC1}   given  \eqn{W}  when   \code{flavor}  is   set  to
### "learning".
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 effIC1,
### The \code{vector}  \code{effIC1} of the  first component of  the efficient
### influence curve, as currently estimated, evaluated at our observations.
 muW,
### The  \code{vector} \code{muW}  of the  conditional expectation  of \eqn{X}
### given \eqn{W}, as currently estimated, evaluated at our observations.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
 verbose=FALSE,
### Prescribes the amount of information  output by the function.  Defaults to
### \code{FALSE}.
 ...
### Additional arguments possibly needed.
) {
  ##seealso<< learnG, learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevTheta
  Z <- (obs[, "X"]-muW) * effIC1;

  obsZ <- cbind(obs, Z=Z);
  verbose && str(verbose, obsZ);

  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnDevMu'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("Z ~", theFormula, "+", theFormula2, sep="")
  } else {
    theFormula <- paste("Z ~", theFormula, sep="")
  } 
  formula <- as.formula(theFormula)
  ## formula <- as.formula(Z~W+I(W^2));


  fit <- glm(formula, data=as.data.frame(obsZ), family=gaussian);
  rm(Z, obsZ);
  if (light) {
    fit <- getLightFit(fit);
  }
  
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;

  return(foo)
### Returns the fitted object.
}


learnDevTheta <- function#Estimation of Cond. Expect. of (Y-thetaXW)^2 Given (X,W) 
### Function   for  the   estimation   of  the   conditional  expectation   of
### \code{(Y-thetaXW)^2}  given  \eqn{(X,W)}  when  \code{flavor}  is  set  to
### "learning".
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 thetaXW,
### The \code{vector} \code{thetaXW} of the conditional expectation of \eqn{Y}
### given \eqn{(X,W)}, as currently estimated, evaluated at our observations.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
verbose=FALSE,
### Prescribes the amount of information  output by the function.  Defaults to
### \code{FALSE}.
 ...
### Additional arguments possibly needed.
) {
  ##seealso<< learnG, learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu
  obsZ <- cbind(obs, Z=thetaXW);
  verbose && str(verbose, obsZ);

  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only 'X' and", paste(varNames, collapse=", "), "in 'learnDevTheta'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula <- paste("I((Y-Z)^2) ~ X+", theFormula, "+ X*(",
                        theFormula, ")", sep="")
  } else {
    theFormula <- paste("I((Y-Z)^2) ~ X+", theFormula, sep="")
  } 
  formula <- as.formula(theFormula)
  ## formula <- as.formula(I((Y-Z)^2)~X*W);
  
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

  return(foo)
### Returns the fitted object.
}

learnG <- function#Estimation of Cond. Prob. of X=x_0 Given W
### Function  for   the  estimation   of  the  conditional   probability  that
### \eqn{X=x_0} (the reference value for \eqn{X}) given \eqn{W}, version based
### on 'glm'.
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the \code{function} \code{tmle.npvi}.
 theX0=0,
### The reference value for \eqn{X}.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
 ...
### Additional arguments possibly needed. 
 ) {
  ##seealso<< learnMuAux, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta
  
  ##
  ## 'glm' version 
  ##
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnG'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("I(X==theX0) ~", theFormula, "+", theFormula2, sep=" ")
  } else {
    theFormula <- paste("I(X==theX0) ~", theFormula, sep="")
  }
  formula <- as.formula(theFormula)
  ## formula <- as.formula(I(X==theX0)~W);

    
  fit <- glm(formula, data=as.data.frame(obs), family="binomial");
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    ## predict(fit, newdata=data.frame(W=W), type="response")
    predict(fit, newdata=data.frame(W), type="response")
  }
  attr(foo, 'fit') <- fit;

  return(foo)
### Returns the fitted object.
}


learnMuAux <- function#Estimation of Cond. Expect. of X Given (X!=x_0, W)
### Function  for the  estimation of  the conditional  expectation  of \eqn{X}
### given \eqn{(X\neq x_0, W)}, version based on 'glm'.
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the  \code{function} \code{tmle.npvi}, where only observations
### with \eqn{X\neq 0} are kept.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
 ...
### Additional arguments possibly needed.
) {
  ##seealso<< learnG, learnTheta, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta
  varNames <- setdiff(colnames(obs), c("X", "Y"))
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only", paste(varNames, collapse=", "), "in 'learnMuAux'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula2 <- paste("I(", varNames, "^2)", collapse=" + ", sep="")
    theFormula <- paste("X ~", theFormula, "+", theFormula2, sep=" ")
  } else {
    theFormula <- paste("X ~", theFormula, sep="")
  }
  formula <- as.formula(theFormula)
  ## formula <- as.formula(X~W+I(W^2));


  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(W) {
    predict(fit, newdata=data.frame(W), type="response");
  }
  attr(foo, 'fit') <- fit;
  
  return(foo)
### Returns the fitted object.
}

learnTheta <- function#Estimation of Cond. Expect. of Y given (X,W)
### Function  for the  estimation of  the conditional  expectation  of \eqn{Y}
### given \eqn{(X, W)}, version based on 'glm'.
(obs,
### The  \code{matrix}  of  observations,  see  for  instance  the  \code{obs}
### argument of the  \code{function} \code{tmle.npvi}, where only observations
### with \eqn{X\neq 0} are kept.
 light=TRUE,
### A  \code{logical},  kept  for   compatibility,  which  should  be  set  to
### \code{TRUE} (its default value). This requires that the result of each fit
### be reduced  in size (for  a faster execution). Currently  implemented only
### for flavor \code{learning}.
 ...
### Additional arguments possibly needed.
) {
  ##seealso<< learnG, learnMuAux, learnCondExpX2givenW, learnCondExpXYgivenW, learnDevG, learnDevMu, learnDevTheta
  varNames <- setdiff(colnames(obs), "Y")
  if (length(varNames)>20) {
    varNames <- varNames[1:20]
    warning(paste("Using only 'X' and", paste(varNames, collapse=", "), "in 'learnTheta'"))
  }
  theFormula <- paste(varNames, collapse=" + ")
  if (length(varNames)<=10) {
    theFormula <- paste("Y ~ X+", theFormula, "+ X*(",
                        theFormula, ")", sep="")
  } else {
    theFormula <- paste("Y ~ X+", theFormula, sep="")
  } 
  formula <- as.formula(theFormula);
  ## formula <- as.formula(Y~X*W);

  
  fit <- glm(formula, data=as.data.frame(obs), family=gaussian);
  if (light) {
    fit <- getLightFit(fit);
  }
  foo <- function(XW) {
    predict(fit, newdata=as.data.frame(XW), type="response");
  }
  attr(foo, 'fit') <- fit;

  return(foo)
### Returns the fitted object.
}


### List of default learning algorithms to use in \code{tmle.npvi} when \code{flavor} is set to "learning".
learningLib <- list(learnCondExpX2givenW=learnCondExpX2givenW,
                    learnCondExpXYgivenW=learnCondExpXYgivenW,
                    learnDevG=learnDevG,
                    learnDevMu=learnDevMu,
                    learnDevTheta=learnDevTheta,
                    learnG=learnG,
                    learnMuAux=learnMuAux,
                    learnTheta=learnTheta)

## run 'makeInstall' with the instruction below commented to make sure the proper Rd files are generated
## then run 'R CMD build' with the instruction below uncommented

rm(learnCondExpX2givenW,    learnCondExpXYgivenW,    learnDevG,    learnDevMu, learnDevTheta, learnG, learnMuAux, learnTheta)



