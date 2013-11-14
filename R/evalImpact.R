evalImpact <- function(res, obs) {
  impactMax <- matrix(NA, ncol=4, nrow=length(res)-1,
                      dimnames=list(paste("step", 1:(length(res)-1), sep=""), c("theta", "mu", "g", "sigma2")))
  impactL1 <- matrix(NA, ncol=4, nrow=length(res)-1,
                     dimnames=list(paste("step", 1:(length(res)-1), sep=""), c("theta", "mu", "g", "sigma2")))
  theta <- sapply(res, function(xx){getTheta(xx)(obs[, c("X", "W")])})
  mu <- sapply(res, function(xx){getMu(xx)(obs[, "W"])})
  g <- sapply(res, function(xx){getG(xx)(obs[, "W"])})
  sigma2 <- sapply(res, getSigma2)
  ##
  dTheta <- apply(theta, 1, diff)
  dMu <- apply(mu, 1, diff)
  dG <- apply(g, 1, diff)
  dSigma2 <- diff(sigma2)
  ## MAX
  fun <- function(xx){max(abs(xx))}
  dt <- apply(dTheta/t(theta[, 1:length(res)-1]), 1, fun)
  dm <- apply(dMu/t(mu[, 1:length(res)-1]), 1, fun)
  dg <- apply(dG/t(g[, 1:length(res)-1]), 1, fun)
  impactMax[, 1] <- dt
  impactMax[, 2] <- dm
  impactMax[, 3] <- dg
  impactMax[, 4] <- dSigma2/sigma2[1:(length(res)-1)]
  ## L1
  fun <- function(xx){mean(abs(xx))}
  dt <- apply(dTheta/t(theta[, 1:length(res)-1]), 1, fun)
  dm <- apply(dMu/t(mu[, 1:length(res)-1]), 1, fun)
  dg <- apply(dG/t(g[, 1:length(res)-1]), 1, fun)
  impactL1[, 1] <- dt
  impactL1[, 2] <- dm
  impactL1[, 3] <- dg
  impactL1[, 4] <- dSigma2/sigma2[1:(length(res)-1)]

  return(list(max=impactMax, L1=impactL1))
}

############################################################################
## HISTORY:
## 2011-01-31
## o Created.
############################################################################

