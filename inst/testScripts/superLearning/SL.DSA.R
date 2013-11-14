SL.DSA <- function (Y, X, newX, family, obsWeights, maxsize = ncol(X),
                    maxorderint = 1, maxsumofpow = 1, Dmove = TRUE, Smove = TRUE,
                    vfold = 5, ...) {
  tryCatch(require(DSA), warning = function(...) {
    stop("you have selected DSA as a library algorithm but do not have the DSA package installed")
  })
  dsaweights <- matrix(obsWeights, nrow = (vfold + 1), ncol = nrow(X),
                       byrow = TRUE)
  fit.DSA <- DSA(Y ~ 1, data = data.frame(Y, X),
                 family = family, maxsize = maxsize, maxorderint = maxorderint,
                 maxsumofpow = maxsumofpow, Dmove = Dmove, Smove = Smove,
                 vfold = vfold, weights = dsaweights)
  pred <- predict(fit.DSA, newdata = newX)
  if (family$family == "binomial") {
    pred <- 1/(1 + exp(-pred))
  }
  fit <- list(object = fit.DSA)
  class(fit) <- "SL.DSA"
  foo <- list(pred = pred, fit = fit)
  return(foo)
}
environment(SL.DSA) <- asNamespace("SuperLearner")

SL.DSA.2 <- function (..., X, maxsize = 2 * ncol(X), maxorderint = 2,
                      maxsumofpow = 2, Smove = FALSE, vfold = 10) {
  SL.DSA(..., X = X, maxsize = maxsize, maxorderint = maxorderint,
         maxsumofpow = maxsumofpow, Smove = Smove, vfold = vfold)
}
environment(SL.DSA.2) <- asNamespace("SuperLearner")


predict.SL.DSA <- function (object, newdata, family, X = NULL, Y = NULL, ...) {
  tryCatch(require(DSA), warning = function(...) {
    stop("you have selected DSA as a library algorithm but do not have the DSA package installed")
  })
  pred <- predict(object = object$object, newdata = newdata)
  if (family$family == "binomial") {
    pred <- 1/(1 + exp(-pred))
  }
  return(pred)
}
environment(predict.SL.DSA) <- asNamespace("SuperLearner")
