this <- npvi
##
##
##
theta <- getTheta(this, tabulate=FALSE)
thetatab <- getTheta(this, tabulate=TRUE)
theta0 <- getTheta0(this, tabulate=FALSE)
theta0tab <- getTheta0(this, tabulate=TRUE)
mu <- getMu(this, tabulate=FALSE)
mutab <- getMu(this, tabulate=TRUE)
g <- getG(this, tabulate=FALSE)
gtab <- getG(this, tabulate=TRUE)
obs <- getObs(this, tabulate=FALSE)
obstab <- getObs(this, tabulate=TRUE)
##
t <- theta(obs[, c("X", "W")])
ttab <- thetatab(obstab[, c("X", "W")])
t0 <- theta0(obs[, "W"])
t0tab <- theta0tab(obstab[, "W"])
##
m <- mu(obs[, "W"])
mtab <- mutab(obstab[, "W"])
##
gg <- g(obs[, "W"])
ggtab <- gtab(obstab[, "W"])
##
##
##

