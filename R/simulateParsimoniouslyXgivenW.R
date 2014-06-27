simulateParsimoniouslyXgivenW <- function(W, obsX, condMeanX, sigma2, parameters, r=3, nMax=10L) {
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Argument 'W':
  W <- Arguments$getNumerics(W);
  if (!is.integer(W)) {
    throw("Parameter 'W' must contain integers entries.")
  }

  ## Argument 'obsX':
  if (sum(obsX==0)!=0) {
    throw("Copy neutral state (obs[,'X']==0) has to be removed");
  }
  
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
  
  getTriangle <- function(a, b, A, B) {
    ## finds points (a1,b1), (a2,b2) and  (a3,b3) from 'cbind(A,B)' such that 
    ## (a,b) belongs to the corresponding triangle
    getBase <- function(ii) {
      ## looks for and returns index 'jj' such that the line going through
      ## (leftA[ii], leftB[ii]) and (rightA[jj], rightB[jj]) is above (a,b),
      ## trying to minimize the distance rightA[jj]-leftA[ii]
      slope <- (b-leftB[ii])/(a-leftA[ii])
      intercept <- (leftB[ii]-leftA[ii]*slope)
      candidates <- which(rightB-slope*rightA-intercept >= 0)
      if (!length(candidates)) {## no candidate
        out <- NA
      } else {
        idx <- which.min(abs(rightA[candidates]-leftA[ii]))
        out <- candidates[idx]
      }
      return(out)
    }
    getThirdVertex <- function(best, opt="L2") {
      ## determines last vertex of triangle
      
      ## first step:
      ## - - - - - - 
      ## looks for and returns indices 'jj' such that the lines going through
      ## (1) left vertex and (a,b), and (2) right vertex and (a,b) are above
      ## points (A[jj], B[jj])
      slopes <- (b-B[best[1:2]])/(a-A[best[1:2]])
      intercepts <- (B[best[1:2]]-A[best[1:2]]*slopes)
      candidates <- which(B-slopes[1]*A-intercepts[1] < 0 &
                          B-slopes[2]*A-intercepts[2] < 0)
      
      ## second step:
      ## - - - - - - -
      pow <- switch(opt, L2=2, L1=1)
      dist <- abs(a-A[candidates])^pow+abs(b-B[candidates])^pow
      best3 <- which.min(dist)
      best3 <- which(A==A[candidates][best3])[1]
      return(best3)
    }
    
    wLeftA <- which(A<a)
    leftA <- A[wLeftA]
    leftB <- B[wLeftA]
    wRightA <- which(A>=a)
    rightA <- A[wRightA]
    rightB <- B[wRightA]
    ## FIXME: lines below
    ##   bases <- cbind(b1=wLeftA, b2=sapply(wLeftA, getBase))
    ## proposal:
    bases <- cbind(b1=seq(along=wLeftA), b2=sapply(seq(along=wLeftA), getBase))
    sLeftA <- sort(leftA, index.return=TRUE)$ix
    sRightA <- sort(rightA, index.return=TRUE)$ix
    crit <- cbind(sLeftA[bases[,1]], sRightA[bases[,2]])
    Best <- bases[which.max(apply(crit, 1, FUN=function(row){
      min(row[1]/length(leftA), 1-row[2]/length(rightA)) ## go as far away from tails as possible
    })),] 

    best <- c(NA, NA, NA)
    best[1] <- which(A==leftA[Best[1]])[1]
    best[2] <- which(A==rightA[Best[2]])[1]
    best[3] <- getThirdVertex(best)
    if (is.na(best[3])) {
      best[3] <- best[1] ## anything would do because given probability 0
    }
    
    return(best)
  }

  getProbs <- function(a, b, A0, B0) {
    if (!is.na(A0[3])) {
      S <- rbind(A0[1:2]-A0[3], B0[1:2]-B0[3])
      probs <- try(solve(S) %*% c(a-A0[3], b-B0[3]))
      if (class(probs)=="try-error") {
        ## browser()
      }
    } else {
      theSum <- sqrt( diff(A0[1:2])^2 + diff(BO[1:2])^2 )
      theDiff <- ( (A0[1]-a)^2 + (BO[1]-b)^2  -
                  (A0[2]-a)^2 + (BO[2]-b)^2 )/theSum
      probs <- 0.5*c(theSum+theDiff, theSum-theDiff)
    }
    out <- c(probs, 1-sum(probs))
    return(out)
  }

  getSimulationScheme <- function(labelW, m1, m2, X) {
    getSimSch <- function(idx) {
      triangle <- getTriangle(m1[idx][1], m2[idx][1], X, X^2)
      if (length(X[triangle])==0) {
        ## browser()
      }
      probs <- getProbs(m1[idx][1], m2[idx][1], X[triangle], X[triangle]^2)
      out <- c(triangle, probs)
      names(out) <- c("X1", "X2", "X3", "p1", "p2", "p3")
      return(out)
    }
    out <- tapply(1:length(labelW), labelW, getSimSch)
    return(out)
  }

  drawFromSimulationScheme <- function(xx, simSch, V) {
    simulationScheme <- simSch[[as.character(labelW[xx[1]])]]
    Xs <- simulationScheme[1:3]
    ps <- simulationScheme[4:6]
    out <- Xs[findInterval(V[xx], cumsum(ps))+1]
    return(out)
  }

  
  phi <- function(x, lambda, xmin=min(obsX), xmax=max(obsX)) {
    lambda*x^2 + (1-lambda)*(x*(xmax+xmin)-xmin*xmax)
  }
  term1 <- (min(obsX)+max(obsX))*mean(parameters$muWB) - mean(1-parameters$gWB)*min(obsX)*max(obsX)
  term2 <- mean(parameters$muWB^2/(1-parameters$gWB))
  lambda <- (term1-sigma2)/(term1-term2)
  if (lambda>1 | lambda<0) {## cannot happen in theory, but may occur due
                            ## to approximations (see the means above)
    lambda <- ifelse(lambda>1, .99, .01)
    cat("Using a slightly modified value for parameter 'lambda' in 'simulateParsimoniouslyXgivenW'...\n")
  }
  condMeanX2 <- phi(condMeanX, lambda)

  obsXq <- unique(quantile(obsX, type=1, probs=seq(0, 1, length=nMax)))
  if (length(setdiff(obsXq, obsX))) {
    throw("This should never happen with type 1 quantiles!")
  }
  tests <- testIfInConvexHull(condMeanX, condMeanX2, obsXq, obsXq^2)
  
  if (FALSE) {
    dev.new()
    xlim <- range(obsXq, condMeanX)
    ylim <- range(obsXq^2, condMeanX2)
    o <- order(obsXq)
    plot(obsXq[o], obsXq[o]^2, xlim=xlim, ylim=ylim, t='l')
    points(condMeanX, condMeanX2, col=2)
  }

  if (!all(tests)) {## if parsimonious method fails (should seldom happen...)
    ## throw("Parsimonious conditional simulation of X given W failed...\n")
    warning("Parsimonious conditional simulation of X given W under a slightly distorted version of the distribution.\nTry a larger 'nMax'...") 
  } else {
    labelW <- identifyUniqueEntries(W)
    simulationSchemes <- getSimulationScheme(labelW, condMeanX, condMeanX2, obsXq)
    V <- runif(length(labelW))
    theXs <- tapply(1:length(labelW), labelW, drawFromSimulationScheme,
                    simSch=simulationSchemes, V=V)
    simulatedXs <- rep(NA, length(labelW))
    for (lab in unique(labelW)) {
      simulatedXs[which(labelW==lab)] <- theXs[[as.character(lab)]]
    }
    out <- simulatedXs
  }
  
  return(out)
}

############################################################################
## HISTORY:
## 2014-02-07
## o Created.
############################################################################

