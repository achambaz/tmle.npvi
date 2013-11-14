require(lodplot) || install.package("lodplot")
data(chrom.bands, "chrom.bands", package = "lodplot")
chrom <- 18
chromdata <- subset(chrom.bands, chrom.bands$chr == chrom)
if (nrow(chromdata) == 0) {
  stop("invalid chromosome:", chrom)
}
lc <- nchar(chromdata$band)
sel <- !(substr(chromdata$band, lc, lc) %in% letters)
chromdata <- chromdata[sel, ]


plot(sin, xlim=c(0, 130))
paint.chromosome(18)

dat <- halfsibscan()
dat18 <- dat[dat$chr==18,]
  
plot.scan(dat18, col="red", with.x=TRUE)
