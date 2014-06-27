##- - - - - 
## learning
##- - - - - 

load("learning.allChromosomes.1-3000.RData")
learning <- list(descr=descr,
                 TMLE=TMLE)
load("learning.allChromosomes.3001-6000.RData")
learning$TMLE <- c(learning$TMLE, TMLE)
load("learning.allChromosomes.6001-9000.RData")
learning$TMLE <- c(learning$TMLE, TMLE)
load("learning.allChromosomes.9001-11314.RData")
learning$TMLE <- c(learning$TMLE, TMLE)

##- - - - - - - - 
## super-learning
##- - - - - - - -

load("superLearning.allChromosomes.1-3000.RData")
superLearning <- list(descr=descr,
                 TMLE=TMLE)
load("superLearning.allChromosomes.3001-6000.RData")
superLearning$TMLE <- c(superLearning$TMLE, TMLE)
load("superLearning.allChromosomes.6001-9000.RData")
superLearning$TMLE <- c(superLearning$TMLE, TMLE)
load("superLearning.allChromosomes.9000-11314.RData")
superLearning$TMLE <- c(superLearning$TMLE, TMLE[-1])

## - - - - - - - - - 
## removing failures
## - - - - - - - - - 

##
## failures due to 'NA' in 'obs'
##

idxL.NA <- sapply(learning$TMLE,
                  function(xx){length(xx)==1 && length(grep("NA", xx))})
idxSL.NA <- sapply(superLearning$TMLE,
                   function(xx){length(xx)==1 && length(grep("NA", xx))})
learning$TMLE <- learning$TMLE[!idxL.NA]
superLearning$TMLE <- superLearning$TMLE[!idxSL.NA]

##
## other failures
##

## <simpleError in findInterval(V[xx], cumsum(ps)): 'vec' contient des valeurs manquantes (NAs)> 
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep("NA", as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep("NA", as.character(xx)))})
learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##  Parsimonious conditional simulation of X given W failed...
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep("Parsimonious", as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep("Parsimonious", as.character(xx)))})
learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##  Parameter 'sigma2' must be positive...
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep("must be positive", as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep("must be positive", as.character(xx)))})
learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

stop()

mean(idxL) # 4%
mean(idxSL) # 1% 

tmle <- learning$TMLE[!idxL]
learning$TMLE <- tmle
tmle <- superLearning$TMLE[!idxSL] 
superLearning$TMLE <- tmle

