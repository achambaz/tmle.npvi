## -----------------------------------------------------------------------------
## setting the different libraries to use when 'flavor' is "h2oEnsembleLearning"
## -----------------------------------------------------------------------------

library(h2oEnsemble)
library(SuperLearner)

EL.library <- c("h2o.glm.wrapper", "h2o.glm.wrapper.00", "h2o.glm.wrapper.05", "h2o.glm.wrapper.10", 
                "h2o.randomForest.wrapper", "h2o.randomForest.1000x100", 
                "h2o.deeplearning.wrapper", "h2o.deeplearning.Rectifier", "h2o.deeplearning.Tanh") 

learnTheta.library <- c(EL.library);

learnG.library <- c(EL.library);

learnMuAux.library <- c(EL.library);

learnDevG.library <- c(EL.library);

learnDevMu.library <- c(EL.library);

learnDevTheta.library <- c(EL.library);

learnCondExpXYgivenW.library <- c(EL.library);

learnCondExpX2givenW.library <- c(EL.library);


### List of default libraries of algorithms to use in \code{tmle.npvi} when \code{flavor} is set to "h2oEnsembleLearning".
h2oEnsembleLearningLib <- list(learnCondExpX2givenW=learnCondExpX2givenW.library,
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
   EL.library);
