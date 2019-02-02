# MFM
R package for simultaneous fine-mapping of genetic associations for several diseases, in a
Bayesian framework that borrows information between the diseases.

## Installation

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

To install, you may first need some R package dependencies, the packages snpStats, GUESSFM, methods, BMA, Rcpp, RcpArmadillo,
parallel, data.table, gtools, fields. E.g., if you don't have gtools, from inside R, do

```R
install.packages("gtools") 
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
``` R 
# R < 3.5:
install.packages("mixOmics")
```

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

