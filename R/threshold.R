threshold <- function(xx, min=0, max=1-min) {
  min <- Arguments$getNumeric(min);
  max <- Arguments$getNumeric(max);
  xx <- Arguments$getNumerics(xx);
  pmin(max, pmax(min, xx))
}

## ##########################################################################
## HISTORY:
## 2011-02-08
## o Changed arguments so that 'threshold' can be used to threshold 'mu'
##   and 'theta' too.
## 2010-12-31
## o Created.
## ##########################################################################

