library(biomaRt);

## assumes that 'geneNames' is defined in global env
## hack:
## load("results/TCGA,OV,expCnMeth,chr18,1-130,2011-09-23.rda")
## geneNames <- row.names(resMat)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annDat <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "band"),
                                filters="hgnc_symbol",
                                values=geneNames,
                                mart=ensembl,
                                uniqueRows=TRUE)
colnames(annDat) <- c("gene", "chromosome", "start", "end", "band");

## remove funky chromosomes
chrs <- annDat$chromosome;
chrs[chrs=="X"] <- "23";
chrs[chrs=="Y"] <- "24";
chrs <- as.numeric(chrs);
idxs <- which(!is.na(chrs));

ad <- annDat[idxs, ];

mm <- match(geneNames, ad$gene)
print(sum(is.na(mm)))

annDat <- ad[mm,]
library(R.utils)
## saveObject(annDat, file="extdata/TCGA,OV,expCnMeth,chr18,annotation,2011-09-29,.xdr")
