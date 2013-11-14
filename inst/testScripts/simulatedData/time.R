## Example of 1+7 TMLE iterations on 200 observations
y <- c(9.82092, 21.41405, 41.4151, 1.111983*60, 1.589712*60, 2.204823*60, 2.856906*60, 3.601251*60)
x <- seq(along=y)-1

## polynomial regression
z <- x[-1]  ## first point is not TMLE
fit <- lm(y[-1]~poly(z, degree=2, raw=TRUE)) 
print(fit)

rhs <- paste(round(fit$coefficients, 1), c("1", "x", "x^2"), collapse="+");
ttl <- paste("y", rhs, sep="=")
xx <- 1:30
yy <- predict(fit, newdata=data.frame(z=xx))

library(R.utils)
pdf("time.pdf", width=20/cm(1), height=20/cm(1))
par(mar=c(5, 5, 1.5, 1), cex.lab=2, cex.axis=2)
plot(xx, yy/60, t='l', xlab="TMLE iteration index", ylab="Time (minutes)")
stext("Time per iteration", side=3, pos=0)
stext(ttl, side=3, pos=1)
points(x, y/60, pch=19)
dev.off()

ttl <- paste("dy/dx", rhs, sep="=")
pdf("cumTime.pdf", width=20/cm(1), height=20/cm(1))
par(mar=c(5, 5, 1.5, 1), cex.lab=2, cex.axis=2)
plot(xx, cumsum(yy)/60, t='l', xlab="TMLE iteration index", ylab="Time (minutes)")
stext("Total time at each iteration", side=3, pos=0)
stext(ttl, side=3, pos=1)
points(x, cumsum(y)/60, pch=19)
dev.off()
