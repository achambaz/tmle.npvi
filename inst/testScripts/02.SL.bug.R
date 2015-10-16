library(tmle.npvi)
example(tmle.npvi)
lib <- lapply(superLearningLib, function(xx){rev(rev(xx)[-1])})
npvi <- tmle.npvi(obs, f=identity, flavor="superLearning", B=5e4, nMax=10, lib=lib)
