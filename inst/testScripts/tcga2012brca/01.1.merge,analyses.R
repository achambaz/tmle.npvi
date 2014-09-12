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

chr <- sapply(names(learning$TMLE),
              function(ll){unlist(strsplit(ll, split=","))[1]})
chr <- sapply(chr, function(ll){unlist(strsplit(ll, split="chr"))[2]})
chr <- as.integer(chr)
posRel <- sapply(names(learning$TMLE),
                 function(ll){unlist(strsplit(ll, split=","))[2]})
posRel <- as.integer(posRel)*1e-3
limChr <- tapply(posRel, chr, max)
cumLimChr <- c(0, cumsum(limChr))

##- - - - - - - - 
## super-learning
##- - - - - - - -

if (FALSE) {
  load("superLearning.allChromosomes.1-3000.RData")
  superLearning <- list(descr=descr,
                        TMLE=TMLE)
  load("superLearning.allChromosomes.3001-6000.RData")
  superLearning$TMLE <- c(superLearning$TMLE, TMLE)
  load("superLearning.allChromosomes.6001-9000.RData")
  superLearning$TMLE <- c(superLearning$TMLE, TMLE)
  load("superLearning.allChromosomes.9001-11314.RData")
  superLearning$TMLE <- c(superLearning$TMLE, TMLE)
  
  superLearning0 <- superLearning
} else {
  load("superLearning.wholeGenome.9001-11314.RData")
  superLearning <- list(descr=descr,
                        TMLE=TMLE)

  superLearning0 <- superLearning
}


## - - - - - - - - - 
## REMOVING FAILURES
## - - - - - - - - - 

## - - - - - - - - - - - - - - - -
## failures due to 'NA' in 'obs'
## - - - - - - - - - - - - - - - -

msg <- "contains at least one 'NA'"
idxL <- sapply(learning$TMLE,
                  function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                  function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 0
sum(idxSL) # 0
learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

## - - - - - - - -
## other failures
## - - - - - - - -

##
## 'simpleError' 
##

msg <- "simpleError"
idxL <- which(sapply(learning$TMLE, inherits, msg))
idxSL <- which(sapply(superLearning$TMLE, inherits, msg))
length(idxL) # 34
length(idxSL) # 2

errL <- sapply(learning$TMLE[idxL], as.character) 
errSL <- sapply(superLearning$TMLE[idxSL], as.character) 

##
## compare to known error messages
##

knownMsgs <- c("Error in findInterval",
               "impossible",
               "Parsimonious conditional simulation of X given W failed",
               "sLeftA",
               "Range of argument",
               "Error in if \\(ps",
               "must be positive")
msgsL <- lapply(knownMsgs, grep, errL)
names(msgsL) <- knownMsgs
str(msgsL)

if (sum(sapply(msgsL, length)) < length(errL)) {
  stop("'learning': some messages are not in 'knownMsgs'")
}
str(msgsL)

msgsSL <- lapply(knownMsgs, grep, errSL)
if (sum(sapply(msgsSL, length)) < length(errSL)) {
  stop("'superLearning': some messages are not in 'knownMsgs'")
}
head(errSL[-unlist(msgsSL)])

##
##  impossible... 
##

msg <- "impossible"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 6 
sum(idxSL) # 0 

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))  # "chr11,002466,KCNQ1"
print(names(superLearning$TMLE[min(which(idxSL))]))  # NA

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]


##
##  Error in if (ps[3]==0)
##

msg <- "Error in if \\("
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 6 
sum(idxSL) # 1 

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))  # "chr11,070245,CTTN"
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr7,056019,GBAS"

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]


##
##  Parameter 'sigma2' must be positive...
##

msg <- "must be positive"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 2 (a weird error!)
sum(idxSL) # 0 

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))  # "chr1,042628,GUCA2A"
print(names(superLearning$TMLE[min(which(idxSL))]))  # NA

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]


##
## 'Exception: Range of argument 'weightsW' is out of range 
##

msg <- "Range of argument 'weightsW'"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 6
sum(idxSL) # 0
learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##
## <simpleError in findInterval(V[xx], cumsum(ps)): 'vec' contient des valeurs manquantes (NAs)>
##

msg <- "Error in findInterval"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})

sum(idxL)   # 0
sum(idxSL)  # 0

## get names of an instance of guilty genes

print(names(learning$TMLE[min(which(idxL))]))        # NA
print(names(superLearning$TMLE[min(which(idxSL))]))  # NA

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##
##  Parsimonious conditional simulation of X given W failed...
##

msg <- "Parsimonious conditional simulation of X given W failed"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==1 && length(grep(msg, as.character(xx)))})
sum(idxL)  # 0
sum(idxSL) # 0

## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))        # NA
print(names(superLearning$TMLE[min(which(idxSL))]))  # NA

learning$TMLE <- learning$TMLE[!idxL]
superLearning$TMLE <- superLearning$TMLE[!idxSL]

##
## <simpleError in sLeftA[bases[, 1]]: type 'list' d'indice incorrect>
##

msg <- "Error in sLeftA\\[bases"
idxL <- sapply(learning$TMLE,
               function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
               ## inherits, "simpleError")
idxSL <- sapply(superLearning$TMLE,
                function(xx){length(xx)==2 && length(grep(msg, as.character(xx)))})
                ## inherits, "simpleError")

sum(idxL)   # 14
sum(idxSL)  # 1


## get names of an instance of guilty genes
print(names(learning$TMLE[min(which(idxL))]))        # "chr1,023833,E2F2"
print(names(superLearning$TMLE[min(which(idxSL))]))  # "chr4,004861,MSX1"

str(learning$TMLE[[min(which(idxL))]])
str(superLearning$TMLE[[min(which(idxSL))]])

learning$TMLE <- learning$TMLE[!idxL]
tmle <- superLearning$TMLE[!idxSL]
tmle <- tmle[!sapply(tmle, is.null)]
superLearning$TMLE <- tmle

saveObject(cumLimChr, file="cumLimChr.xdr")
saveObject(learning, file="learningAllChromosomes.xdr")
saveObject(superLearning, file="superLearningAllChromosomes.xdr")
