set.seed(1)
##
## 'identifyUniqueEntries'
##
W <- sample(round(runif(3), 2), 12, replace=TRUE)
labelW <- identifyUniqueEntries(W)
##
## 'getTriangle' and 'getProbs'
##
a <- 1/2
b <- 1/3                                        
A <- runif(1e3)
B <- A^2
##
plot(A, B)
test <- in.chull(a, b, A, B)
if (test) {## if (a,b) in convex hull of A and B
  idx <- getTriangle(a, b, A, B)
  points(A[idx], B[idx], bg=3, pch=21)
  lines(A[rep(idx, 2)], B[rep(idx, 2)], col=3)
  points(a, b, bg=2, pch=21)
  ##
  probs <- getProbs(a, b, A[idx], B[idx])
  V <- A[idx][findInterval(runif(1e4), cumsum(probs))+1]
  table(V)
  print(c(mean(V), mean(V^2)))
} else{## cannot work
  cat("(a,b) not in convex hull of A and B...\n")
}
##
## 'getSimulationScheme'
##
W <- sample(1:2, length(A), replace=TRUE)
labelW <- identifyUniqueEntries(W)
obs <- matrix(A, ncol=1)
colnames(obs) <- "X"
m1 <- c(1/2, 1/3)[labelW]
m2 <- c(1/3, 1/4)[labelW]
#### x11()
plot(A, B)
points(m1, m2, bg=c(2, 3)[labelW], pch=21)
test <- testIfInConvexHull(m1, m2, A, B)
if (all(test)) {
  X <- obs[, "X"]
  X <- X[X!=0]
  simSch <- getSimulationScheme(labelW, m1, m2, X)
  lines(A[rep(simSch[[1]][1:3], 2)],
        B[rep(simSch[[1]][1:3], 2)], col=2)
  points(A[simSch[[1]][1:3]],
         B[simSch[[1]][1:3]], bg=2, pch=21)
  lines(A[rep(simSch[[2]][1:3], 2)],
        B[rep(simSch[[2]][1:3], 2)], col=3)
  points(A[simSch[[2]][1:3]],
         B[simSch[[2]][1:3]], bg=3, pch=21)
  ##
  V <- runif(length(labelW))
  theXs <- tapply(1:length(labelW), labelW, drawFromSimulationScheme,
                  simSch=simSch, V=V, XX=X)
  simulatedXs <- rep(NA, length(labelW))
  for (lab in unique(labelW)) {
    simulatedXs[which(labelW==lab)] <- theXs[[as.character(lab)]]
  }
}
tapply(simulatedXs, labelW, function(xx){c(mean(xx), mean(xx^2))})
