library("R.utils");
log <- Arguments$getVerbose(-8, timestamp=TRUE);

path <- "inst/extdata";
path <- Arguments$getReadablePath(path);
filenames <- list.files(path, pattern="*.txt");

nf <- length(filenames);
if (nf==0) {
  throw("No data file found");
}
if ((nf>1) && require("R.menu")) {
  filename <- textMenu(filenames, value=TRUE);
} else {
  filename <- filenames[1];
}

pathname <- file.path(path, filename);
obs <- read.table(pathname, header=TRUE);

############################################################################
# HISTORY:
# 2010-07-18
# o Created.
############################################################################
