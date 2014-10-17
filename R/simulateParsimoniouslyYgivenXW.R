simulateParsimoniouslyYgivenXW <- function(T, Y) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'T':
  T <- Arguments$getNumerics(T);
  
  ## Argument 'Y':
  Y <- Arguments$getNumerics(Y);
  

  getSimulationScheme <- function(T, Y) {
    ## Yinf <- sapply(T, function(tt) {
    ##   max(Y[which(Y<tt)])
    ## })
    ## Ysup <- sapply(T, function(tt) {
    ##   min(Y[which(Y>=tt)])
    ## })
    sortedY <- sort(Y, index.return=TRUE)
    index <- findInterval(T, sortedY$x)
    Yinf.value <- sortedY$x[index]
    Ysup.value <- sortedY$x[index+1]
    Yinf.index <- sortedY$ix[index]
    Ysup.index <- sortedY$ix[index+1]
    prob <- (Ysup.value - T)/(Ysup.value-Yinf.value)
    out <- cbind(Yinf.index, Ysup.index, prob, 1-prob)
    colnames(out) <- c("Y1", "Y2", "p1", "p2")
    return(out)
  }
  
  test <- (min(Y)<=min(T) & max(T)<=max(Y))
  if (!test) {## if parsimonious method fails
    throw("Parsimonious conditional simulation of Y given (X, W) only works when 'theta' takes its values in 'range(Y)'\nPlease use options 'thetamin' and 'thetamax' in the construction of the 'NPVI' object")
  } else {
    simulationSchemes <- getSimulationScheme(T, Y)
    V <- (runif(length(T)) >= simulationSchemes[, "p1"]) + 1
    out <- simulationSchemes[cbind(1:length(T), V)]
  }
  
  return(out)
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

