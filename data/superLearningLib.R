requireNamespace("SuperLearner")
## -----------------------------------------------------------------------
## setting the different libraries to use when 'flavor' is "superLearning"
## -----------------------------------------------------------------------

## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA")
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm")
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam")

requireNamespace("SuperLearner")
## library(e1071)
## library(DSA)
## library(glmnet)

learnTheta.library <- c("SL.glm.theta", "SL.polymars", SL.library);

learnG.library <- c("SL.glm.g", SL.library);

learnMuAux.library <- c(SL.library);

learnDevG.library <- c(SL.library)

learnDevMu.library <- c(SL.library)

learnDevTheta.library <- c(SL.library)

learnCondExpXYgivenW.library <- c("SL.glm.condExpXYgivenW", SL.library)

learnCondExpX2givenW.library <- c("SL.glm.condExpX2givenW", SL.library)


### List of default libraries of algorithms to use in \code{tmle.npvi} when \code{flavor} is set to "superLearning".
superLearningLib <- list(learnCondExpX2givenW=learnCondExpX2givenW.library,
                         learnCondExpXYgivenW=learnCondExpXYgivenW.library,
                         learnDevG=learnDevG.library,
                         learnDevMu=learnDevMu.library,
                         learnDevTheta=learnDevTheta.library,
                         learnG=learnG.library,
                         learnMuAux=learnMuAux.library,
                         learnTheta=learnTheta.library)


rm(learnTheta.library,
   learnG.library,
   learnMuAux.library,
   learnDevTheta.library,
   learnDevG.library,
   learnDevMu.library,
   learnCondExpXYgivenW.library,
   learnCondExpX2givenW.library,
   SL.library)
