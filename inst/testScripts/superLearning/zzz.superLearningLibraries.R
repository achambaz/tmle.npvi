SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA")##antoine##, "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");

learnTheta <- c("SL.glm.theta", "SL.polymars", SL.library);
## learnTheta0 <- c("SL.glm.theta", "SL.polymars", SL.library);  ## not tested 
learnG <- c("SL.glm.g", SL.library);
## learnMuAux <- c("SL.glm.muAux", SL.library)
learnMuAux <- c(SL.library);

## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm", "SL.DSA", "SL.glmnet");
SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.svm")##antoine##, "SL.glmnet");
## SL.library <- c("SL.glm", "SL.randomForest", "SL.gam", "SL.glmnet");

learnDevG <- c(SL.library)
learnDevMu <- c(SL.library)
learnDevTheta <- c(SL.library)

learnCondExpXYgivenW <- c("SL.glm.condExpXYgivenW", SL.library)
learnCondExpX2givenW <- c("SL.glm.condExpX2givenW", SL.library)

SL.library <- unique(c(learnTheta, learnG, learnMuAux,
                       learnDevTheta, learnDevG, learnDevMu,
                       learnCondExpXYgivenW, learnCondExpX2givenW))
