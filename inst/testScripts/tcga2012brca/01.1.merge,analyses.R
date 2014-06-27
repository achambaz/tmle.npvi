##- - - - - 
## learning
##- - - - - 
## wd: 'clusterData'
load("learning.allChromosomes.1-3000.RData")
learning <- list(descr=descr,
                 TMLE=TMLE)
load("learning.allChromosomes.3001-6000.RData")
learning$TMLE <- c(learning$TMLE, TMLE)
load("learning.allChromosomes.6001-9000.RData")
learning$TMLE <- c(learning$TMLE, TMLE)
load("learning.allChromosomes.9001-11314.RData")
learning$TMLE <- c(learning$TMLE, TMLE)
learning0 <- learning

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

superLearning0 <- superLearning

## - - - - - - - - - 
## removing failures
## - - - - - - - - - 

##
## failures due to 'NA' in 'obs'
##

idxL.NA <- sapply(learning$TMLE,
                  function(xx){length(xx)==1 && length(grep("contains at least one 'NA'", as.character(xx)))})
idxSL.NA <- sapply(superLearning$TMLE,
                  function(xx){length(xx)==1 && length(grep("contains at least one 'NA'", as.character(xx)))})
mean(idxL.NA)  # 9%
mean(idxSL.NA) # 5%
learning$TMLE <- learning$TMLE[!idxL.NA]
superLearning$TMLE <- superLearning$TMLE[!idxSL.NA]

##
## other failures
##

## 1- 'simpleError'
idxL <- sapply(learning$TMLE, inherits, "simpleError")
idxSL <- sapply(superLearning$TMLE, inherits, "simpleError")
sum(idxL)
sum(idxSL)

errL <- sapply(learning$TMLE[idxL], as.character) 
errSL <- sapply(learning$TMLE[idxSL], as.character) 

## compare to known error messages
knownMsgs <- c("Error in findInterval", "impossible", "Parsimonious conditional simulation of X given W failed", "sLeftA")
msgsL <- lapply(knownMsgs, grep, errL)
names(msgsL) <- knownMsgs
str(msgsL)


if (sum(sapply(msgsL, length)) < length(errL)) {
  stop("Some messages are not in 'knownMsgs'")
}
str(msgsL)

msgsSL <- lapply(knownMsgs, grep, errSL)
if (sum(sapply(msgsSL, length)) < length(errSL)) {
  stop("Some messages are not in 'knownMsgs'")
}
head(errSL[-unlist(msgsSL)])




## <simpleError in findInterval(V[xx], cumsum(ps)): 'vec' contient des valeurs manquantes (NAs)> 
msg <- "Error in findInterval"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})

sum(idxL)   # 39
sum(idxSL)  # 4

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))        # "chr10,113910,GPAM"
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr11,007992,EIF3F"

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]


##  Parsimonious conditional simulation of X given W failed...

msg <- "Parsimonious conditional simulation of X given W failed"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
sum(idxL)  # ~ 476 %
sum(idxSL) # ~ 104 %

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))        # "chr10,001095,WDR37"
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr10,017632,PTPLA"

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##  Parameter 'sigma2' must be positive...
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep("must be positive", as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep("must be positive", as.character(xx)))})
sum(idxL)  # 0
sum(idxSL) # 1 (a weird error!)

## get names of an instance of guilty genes
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr10,064927,JMJD1C"

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

## <simpleError in findInterval(V[xx], cumsum(ps)): 'vec' contient des valeurs manquantes (NAs)> 
msg <- "Error in findInterval"
idxL <- sapply(learning$TMLE, inherits, "simpleError")
idxSL <- sapply(superLearning$TMLE, inherits, "simpleError")

sum(idxL)   # 15
sum(idxSL)  # 961

learning$TMLE[idxL]


## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))        # "chr10,113910,GPAM"
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr11,007992,EIF3F"



str(learning$TMLE[[913]])

stop()

tmle <- learning$TMLE[!idxL]
learning$TMLE <- tmle
tmle <- superLearning$TMLE[!idxSL] 
superLearning$TMLE <- tmle

