Package: tmle.npvi
===================

Version: 0.11.5 [2017-02-17]

o Now importing from "SuperLearner" + corresponding updates to NAMESPACE

Version: 0.11.4 [2017-02-17]

o Fixes in NAMESPACE (now importing 'R.utils::cat' and exporting 'as.character')

Version: 0.11.3 [2017-02-16]
o Most of the documentation is now generated from roxygen2 and not inlinedocs
o Added man files to version control

Version: 0.11.0 [2015-10-09]
o The package now handles weights and cluster identification.
o New flavor "h2oEnsembleLearning" added for super learning.

Version: 0.10.0  [2015-05-22]
o Added CITATION file.
o Added reference to Bioinformatics Application Note.
o Enhanced the parsimonious simulation  of 'X' given 'W' thanks
to a discussion with Emily Chang (UCSF). Side effect: package does not
need to import 'sgeostat::in.chull' and 'sgeostat::in.polygon' anymore.
o Passes R CMD check.
o The function 'tmle.npvi' now also handles data frames.

Version: 0.9.3 [2015-02-05]
o Shortened package title as per the CRAN policies.

Version: 0.9.2 [2015-02-05]
o Updated NAMESPACE to make use of S3 method registration.
o Dropped support of 'DSA' SuperLearning libraries, as DSA is not a
'mainstream package'.

Version: 0.9.1 [2014-12-12]
o Added 'getPValue' for testing ``\eqn{Psi(P_0)=Phi(P_0)}'' or ``\eqn{Psi(P_0)=0}''.

Version: 0.9.0 [2014-11-14]
o Faster version, which handles much larger data sets:
  - Now using sparse matrices to speed up the computations of the tabulated
versions of the features of interest. Based on new parameter 'nMax' in 'getSimulationScheme'.
  - Faster version of 'simulateParsimiouslyYgivenXW' and 'simulateParsimiouslyXgivenW'
  - Added parameter 'nMax' to speed up 'getSimulationScheme'.
o Package does not DEPEND on 'MASS' and 'sgeostat' anymore, but IMPORTS
functions from them.
o Added tcga12brca data set.
o Updated DESCRIPTION as suggested by Henrik Bengtsson.
o Enhanced scripts for TCGA data analysis.

Version: 0.8.1 [2014-02-08]
o Now using lazy-loading of learning libraries instead of assigning
  objects to the global environment.

Version: 0.8.0 [2014-02-07]
o Passes R CMD check --as-cran

(...)

Version: 0.1.0 [2010-07-13]
o Passes R CMD check.
o Created.

