fasterGetSimulationScheme <- function(labelW, condMeanX, condMeanX2, Xq.value) {
  ## preliminary
  Xq <- Xq.value
  keepOnly <- match(unique(labelW), labelW)
  lab <- labelW[keepOnly]
  condMeanX <- condMeanX[keepOnly]
  condMeanX2 <- condMeanX2[keepOnly]
  ## preparing triangles
  idx <- 1:length(Xq)
  triangles <- as.matrix(expand.grid(idx, idx, idx))
  keep <- apply(triangles, 1, FUN=function(x) all(diff(x)>0))
  triangles <- triangles[keep, ]
  
  ## ordering triangles by distance to tails
  left <- Xq[triangles[, 1]]-min(Xq)
  right <- max(Xq)-Xq[triangles[, 3]]
  oo <- order(pmin(left, right), decreasing=TRUE)
  triangles <- triangles[oo, ]

  ## assigning a triangle to each couple '(condMeanX[ii], condMeanX2[ii])'
  ## and completing the simulation scheme
  trg <- matrix(NA, ncol=length(condMeanX), nrow=3)
  probs <- matrix(NA, nrow=length(condMeanX), ncol=3)
  for (ii in 1:nrow(triangles)) {
    idx <- triangles[ii, ]
    jdx <- is.na(trg[1, ])
    if (length(jdx)==0) {
      break
    }
    test <- in.polygon(condMeanX[jdx], condMeanX2[jdx],
                       Xq[idx], Xq[idx]^2)
    if (any(test)) {
      concerned <- which(jdx)[which(test)]
      trg[, concerned] <- idx
      probs[concerned, ] <- cart2bary(cbind(Xq[idx], Xq[idx]^2),
                                      cbind(condMeanX[concerned], condMeanX2[concerned]))
    }
  }
  out <- cbind(t(trg), probs)
  colnames(out) <- c("X1", "X2", "X3", "p1", "p2", "p3")
  ## out <- lapply(seq_len(nrow(out)), function(ii){out[ii, ]})
  ## names(out) <- as.character(lab)
  ## out <- out[order(lab)]
  rownames(out) <- as.character(lab)

  return(out)
}

simulateParsimoniouslyXgivenW <- function(W, xmin, xmax, Xq, condMeanX, sigma2, parameters, r=3) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'W':
  W <- Arguments$getNumerics(W);
  if (!is.integer(W)) {
    throw("Parameter 'W' must contain integers entries.")
  }

  ## Arguments 'xmin' and 'xmax':
  xmin <- Arguments$getNumerics(xmin)
  xmax <- Arguments$getNumerics(xmax)
  if (xmin>=xmax) {
    throw("Argument 'xmin' must be smaller than argument 'xmax'...")
  }
  
  ## Argument 'Xq':
  Xq.value <- Arguments$getNumerics(Xq$value)
  Xq.index <- Arguments$getNumerics(Xq$index)
  ## ## CAUTION!
  Xq.value <- Xq.value[Xq.value!=0]
  
  ## Argument 'condMeanX':
  condMeanX <- Arguments$getNumerics(condMeanX);
  
  ## Argument 'sigma2':
  sigma2 <- Arguments$getNumeric(sigma2);
  if (sigma2<=0) {
    throw("Parameter 'sigma2' must be positive.")
  }

  if (length(W) != length(condMeanX)) {
    throw("Vectors 'W' and 'condMeanX' must be of same length.")
  }
  

  ## Argument 'parameters'
  ## TODO 
  
  ## Argument 'r':
  r <- Arguments$getInteger(r);
  if (r!=3) {
    throw("Only the case 'r=3' is implemented.");
  }

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## A few useful functions
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  identifyUniqueEntries <- function(W) {
    ## attributes a unique label to each unique entry of W
    sW <- sort(W, index.return=TRUE)
    labelSW <- cumsum(c(0, diff(sW$x)>0))
    labelW <- rep(NA, length(W))
    labelW[sW$ix] <- labelSW
    return(labelW+1)
  }

  testIfInConvexHull <- function(a, b, A, B) {
    ## tests if points with abscissa and ordinate of the form 'a[ii]' and 'b[i]' belong
    ## to convex hull of points with abscissa  and ordinate of the form 'A[jj]' and 'B[jj]'
    ab <- unique(cbind(a, b))
    in.chull(ab[, 1], ab[, 2], A, B)
  }
  
  
  phi <- function(x, lambda, x.min=xmin, x.max=xmax) {
    lambda*x^2 + (1-lambda)*(x*(x.max+x.min)-x.min*x.max)
  }
  term1 <- (xmin+xmax)*mean(parameters$muWB) - mean(1-parameters$gWB)*xmin*xmax
  term2 <- mean(parameters$muWB^2/(1-parameters$gWB))
  lambda <- (term1-sigma2)/(term1-term2)
  if (lambda>1 | lambda<0) {## cannot happen in theory, but may occur due
                            ## to approximations (see the means above)
    lambda <- ifelse(lambda>1, .99, .01)
    cat("Using a slightly modified value for parameter 'lambda' in 'simulateParsimoniouslyXgivenW'...\n")
  }
  condMeanX2 <- phi(condMeanX, lambda)

  tests <- testIfInConvexHull(condMeanX, condMeanX2, Xq.value, Xq.value^2)
  
  if (FALSE) {
    dev.new()
    xlim <- range(Xq.value, condMeanX)
    ylim <- range(Xq.value^2, condMeanX2)
    o <- order(Xq.value)
    plot(Xq.value[o], Xq.value[o]^2, xlim=xlim, ylim=ylim, t='l')
    points(condMeanX, condMeanX2, col=2)
  }

  if (!all(tests)) {## if parsimonious method fails (should seldom happen...)
    ## throw("Parsimonious conditional simulation of X given W failed...\n")
    warning("Parsimonious conditional simulation of X given W under a slightly distorted version of the distribution. You may want to try a larger 'nMax'...") 
  } 
  labelW <- identifyUniqueEntries(W)
  simulationSchemes <-  fasterGetSimulationScheme(labelW, condMeanX, condMeanX2, Xq.value)
  V <- runif(length(labelW))
  idx <- match(labelW, unique(labelW))
  sch <- simulationSchemes[idx, ]
  Xs <- sch[, 1:3]
  ps <- sch[, 4:6]
  q1 <- ps[, 1]
  q2 <- ps[, 2] + q1
  randomIndex <- 1+(V>q1)+(V>q2)
  simulatedXs <- Xs[cbind(1:length(randomIndex), randomIndex)]

  out <- simulatedXs
  return(out)
}



############################################################################
## HISTORY:
## 2014-02-07
## o Created.
## 2014-11-14
## o Substantial speedup by avoiding 'tapply'.
############################################################################

