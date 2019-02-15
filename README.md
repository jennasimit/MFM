# MFM
R package for simultaneous fine-mapping of genetic associations for several diseases, in a
Bayesian framework that borrows information between the diseases.

Website available at: https://jennasimit.github.io/MFM/.

MFMextra has simulation code for evaluating MFM and is available at: https://jennasimit.github.io/MFMextra/.
A MFM demo is available in the vignettes of MFM. In our tests in ran in under 5 seconds.

## System Requirements

MFM could be installed with ease on versions of R > 3.5 and requires additional steps for R < 3.5.
Installation has been tested on R 3.5.1. Installation time is estimated as 2 minutes.

## Installation Guide

### Short version

```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("jennasimit/MFM")
```

### Longer version (if this fails)

MFM is purely an R package and so is platform independent. It depends on the R package GUESSFM, which depends on the 
software GUESS, available at http://www.bgx.org.uk/software/guess.html, but can be more easily installed via the R package 
R2GUESS, which is a dependency of GUESSFM.
GUESSFM is described in the paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005272 and 
GUESS is described in the paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003657.

To install, you may first need some R package dependencies, the packages snpStats, GUESSFM, BMA, Rcpp, RcpArmadillo,
parallel, data.table, gtools, fields, cowplot. E.g., if you don't have some of these packages, from inside R, do

```R
install.packages("BMA") 
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("parallel")
install.packages("data.table")
install.packages("gtools")
install.packages("fields")
install.packages("cowplot")
```
 
Some packages (e.g. snpStats) are from Bioconductor.  For these, you need to do

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("snpStats")
```
R2GUESS, a dependency of GUESSFM, has the dependency mixOmics, which is available from Bioconductor for R versions >=3.5
and at cran in archive form for R < 3.5.
``` R 
# R >= 3.5:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("mixOmics")
```

R < 3.5:

Go to https://cran.r-project.org/src/contrib/Archive/mixOmics/ and download mixOmics_6.3.2.tar.gz. Then, in R:
```R
install.packages("mixOmics", repos = "mixOmics_6.3.2.tar.gz", type="source")
```
where if in a different working directory than where the zip file is contained, a file path may be needed for the repos
argument.

GUESSFM is from github and to install do
```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("chr1swallace/GUESSFM")
```

Then, to install MFM, also from github, do
```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("jennasimit/MFM")
```


Run the following to generate the MFM website:
```
Rscript -e "pkgdown::build_site()"
```

