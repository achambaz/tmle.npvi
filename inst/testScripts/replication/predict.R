## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Check that predict gives identical values before and after 'getLightFit'
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N <- 1e5
x <- rnorm(N)
y <- 3*x + 1 + rnorm(N)
xx <- rnorm(1e3)

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Fit a linear model
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- lm(y~x)
## save(fit, file="fit.rda")

lfit <- getLightFit(fit)
## save(lfit, file="lfit.rda")

pred <- predict(fit, newdata=data.frame(x=xx))
lpred <- predict(lfit, newdata=data.frame(x=xx))
stopifnot(all.equal(pred, lpred))


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Fit a generalized linear model
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- glm(y~x)
## save(fit, file="fit.rda")

lfit <- getLightFit(fit)
## save(lfit, file="lfit.rda")

pred <- predict(fit, newdata=data.frame(x=xx))
lpred <- predict(lfit, newdata=data.frame(x=xx))
stopifnot(all.equal(pred, lpred))


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Fit a 'rpart' model
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit <- rpart(y~x)
## save(fit, file="fit.rda")

lfit <- getLightFit(fit)
## save(lfit, file="lfit.rda")

pred <- predict(fit, newdata=data.frame(x=xx))
lpred <- predict(lfit, newdata=data.frame(x=xx))
stopifnot(all.equal(pred, lpred))
