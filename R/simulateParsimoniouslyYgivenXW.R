simulateParsimoniouslyYgivenXW <- function(T, Y) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'T':
  T <- Arguments$getNumerics(T);
  
  ## Argument 'Y':
  Y <- Arguments$getNumerics(Y);
  
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## A few useful functions
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  identifyUniqueEntries <- function(T) {
    ## attributes a unique label to each unique entry of T
    sT <- sort(T, index.return=TRUE)
    labelST <- cumsum(c(0, diff(sT$x)>0))
    labelT <- rep(NA, length(T))
    labelT[sT$ix] <- labelST
    return(labelT+1)
  }

  getSimulationScheme <- function(labelT, T, Y) {
    getSimSch <- function(idx) {
      t <- T[idx][1]
      wLeftY <- which(Y<=t)
      wRightY <- which(Y>=t)
      leftY <- max(Y[wLeftY])
      rightY <- min(Y[wRightY])
      if (length(wRightY)==0) {
        warning("Parsimonious conditional simulation of Y given (X, W) only works when 'theta' takes its values in 'range(Y)'")
        prob <- 0
        out <- c(which(Y==leftY)[1], t, prob, 1-prob)
      } else if (length(wLeftY)==0) {
        warning("Parsimonious conditional simulation of Y given (X, W) only works when 'theta' takes its values in 'range(Y)'")
        prob <- 1
        out <- c(t, which(Y==rightY)[1], prob, 1-prob)
      } else {
        prob <- (rightY-t)/(rightY-leftY)  ## is NA if 'rightY==leftY'
        out <- c(which(Y==leftY)[1], which(Y==rightY)[1], prob, 1-prob)
      }
      names(out) <- c("Y1", "Y2", "p1", "p2")
      return(out)
    }
    out <- tapply(1:length(labelT), labelT, getSimSch)
    return(out)
  }

  
  drawFromSimulationScheme <- function(xx, simSch, V) {
    simulationScheme <- simSch[[as.character(labelT[xx[1]])]]
    Ys <- as.integer(simulationScheme[1:2])
    if (Ys[1]==Ys[2]) {
      out <- rep(Ys[1], length(xx))
    } else {
      ps <- simulationScheme[3:4]
      out <- Ys[findInterval(V[xx], cumsum(ps))+1]
    }
    return(out)
  }

  test <- (min(Y)<=min(T) & max(T)<=max(Y))
  if (!test) {## if parsimonious method fails
    throw("Parsimonious conditional simulation of Y given (X, W) only works when 'theta' takes its values in 'range(Y)'\nPlease use options 'thetamin' and 'thetamax' in the construction of the 'NPVI' object")
  } else {
    labelT <- identifyUniqueEntries(T)
    simulationSchemes <- getSimulationScheme(labelT, T, Y)
    V <- runif(length(labelT))
    theYs <- tapply(1:length(labelT), labelT, drawFromSimulationScheme,
                    simSch=simulationSchemes, V=V)
    simulatedYs <- rep(NA, length(labelT))
    for (lab in unique(labelT)) {
      simulatedYs[which(labelT==lab)] <- theYs[[as.character(lab)]]
    }
    out <- simulatedYs
  }
  
  return(out)
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

