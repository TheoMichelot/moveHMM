# moveHMM
An R package for animal movement modelling using hidden Markov models.

Get started with the vignette: [Guide to using moveHMM](https://cran.r-project.org/web/packages/moveHMM/vignettes/moveHMM-guide.pdf)

## Installation instructions

### Stable release
The package is available on [CRAN](https://cran.r-project.org/web/packages/moveHMM/index.html). To install it from CRAN,
you can use the following commands:
``` R
# install dependencies
install.packages(c("Rcpp","RcppArmadillo","sp","CircStats"))
# install moveHMM
install.packages("moveHMM")
```

### Install from Github
To install the latest (**unstable**) version of the package from Github:
``` R
install.packages(c("Rcpp","RcppArmadillo","sp","CircStats","devtools"))
library(devtools)
install_github("TheoMichelot/moveHMM", build_vignettes=TRUE)
```
