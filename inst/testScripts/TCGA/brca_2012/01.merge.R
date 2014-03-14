library(R.utils)
path <- "data"
path <- Arguments$getReadablePath(path)

filename <- "BRCA.methylation.27k.450k.466.txt"
pathname <- file.path(path, filename)

mdat <- read.table(pathname, sep="\t", check.names=FALSE)
mpnames <- names(mdat)
mpnames <- substr(mpnames, 1, 15)
str(mpnames)
names(mdat) <- mpnames
methMat <- as.matrix(mdat)

## annotation data (methylation)
apath <- "annotationData/chipTypes/HumanMethylation450"
apath <- Arguments$getReadablePath(apath)
afilename <- "HumanMethylation450_15017482_v.1.2.csv"
apathname <- file.path(apath, afilename)

adat <- read.table(apathname, sep=",", skip=7, quote="\"", nr=485577, header=TRUE, as.is=TRUE)
names(adat)

## focus on probes that are on the microarray
cgNames <- rownames(mdat)  ## names of probes on the microarray
mm <- match(cgNames, adat[["Name"]])
adat2 <- adat[mm, ]
geneNames <- adat2[, "UCSC_RefGene_Name"]
mgnames <- unique(unlist(strsplit(geneNames, ";")))

## painful effort 
x <- strsplit(geneNames, ";")
x <- lapply(x, unique)
xx <- sapply(x, paste, collapse=";")
length(xx)
write.table(xx, file="toto.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)
system("./findGenesIndices.pl > tata.txt")
y <- read.table("tata.txt", sep="\t", header=FALSE, as.is=TRUE)
names(y) <- c("name", "char")
head(y)

dim(y)
length(intersect(y$name, gids))==length(gids)

mg <- match(gids, y$name)

if (FALSE) { ## sanity check
  ids <- as.numeric(unlist(strsplit(y$char[mg[1]], " ")))
  unique(unlist(x[ids]))==gids[1]

  lens <- sapply(y$char, FUN=function(x) length(unlist(strsplit(x, " "))))
  hist(lens)
}

## copy number gene names
cpathname <- "data/GISTIC2/all_data_by_genes.txt"
cdat <- read.table(cpathname, sep="\t", as.is=TRUE, header=TRUE, check.names=FALSE)
cnMat <- as.matrix(cdat[, -c(1:3)])

cpnames <- colnames(cnMat)
## cpnames <- substr(cpnames, 1, 15)
str(cpnames)

cgnames <- cdat[["Gene Symbol"]]
str(cgnames)

## expression gene and patient names
epathname <- "data/BRCA.exp.466.med.txt"
edat <- read.table(epathname, sep="\t", as.is=TRUE, header=TRUE, quote="\"", check.names=FALSE)
egnames <- edat[["NAME"]]
str(egnames)

exprMat <- as.matrix(edat[, -1])
colnames(exprMat) <- substr(colnames(exprMat), 1, 15)
epnames <- colnames(exprMat)
str(epnames)

## gene names
length(intersect(cgnames, mgnames))  ## 13261!
length(intersect(cgnames, egnames))  ## 15466!
length(intersect(intersect(cgnames, mgnames), egnames))  ## 11943

gids <- intersect(intersect(cgnames, mgnames), egnames)

eg <- match(gids, egnames)
cg <- match(gids, cgnames)

## patient names
length(intersect(cpnames, mpnames)) 
length(intersect(cpnames, epnames)) 
length(intersect(intersect(cpnames, mpnames), epnames))

pids <- intersect(intersect(cpnames, mpnames), epnames)

ep <- match(pids, epnames)
cp <- match(pids, cpnames)
mp <- match(pids, mpnames)

## export
gid <- gids[1]

## gene expression
idxE <- eg[match(gid, gids)]
geneExpr <- exprMat[idxE, ep]
str(geneExpr)

## DNA copy number
idxC <- cg[match(gid, gids)]
copyNumber <- cnMat[idxC, cp]
str(copyNumber)

stopifnot(identical(names(copyNumber), names(geneExpr)))

## methylation
idxM <- mg[match(gid, gids)]
idxsM <- as.numeric(unlist(strsplit(y[idxM, 2], " ")))
methyl <- methMat[idxsM, mp, drop=FALSE]
dim(methyl)
stopifnot(identical(names(copyNumber), colnames(methyl)))

rownames(methyl) <- paste("W", 1:nrow(methyl), sep="")

obs <- cbind(Y=geneExpr, X=copyNumber, W=t(methyl))
str(obs)

pairs(obs)


saveObject(obs, file="obs,ATAD3A.xdr")

