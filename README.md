
# moveHMM [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/moveHMM)](https://cran.r-project.org/package=moveHMM) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) [![Downloads](http://cranlogs.r-pkg.org/badges/moveHMM)](https://cran.r-project.org/package=moveHMM)

An R package for modelling animal movement with hidden Markov models.

Get started with the vignette: [Guide to using moveHMM](https://CRAN.R-project.org/package=moveHMM/vignettes/moveHMM-guide.pdf)

## Installation instructions

### Stable release
The package is available on [CRAN](https://CRAN.R-project.org/package=moveHMM). To install it from CRAN, you can use the following commands:
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

## References
Michelot, T., Langrock, R., Patterson, T.A. (2016). [moveHMM: An R package for analysing animal movement data using hidden Markov models](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12578/abstract). *Methods in Ecology and Evolution*. 7(11), 1308-1315.

Langrock, R., King, R., Matthiopoulos, J., Thomas, L., Fortin, D., & Morales, J. M. (2012). [Flexible and practical modeling of animal telemetry data: hidden Markov models and extensions](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/11-2241.1). *Ecology*, 93(11), 2336-2342.

Patterson, T. A., Basson, M., Bravington, M. V., & Gunn, J. S. (2009). [Classifying movement behaviour in relation to environmental conditions using hidden Markov models](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2656.2009.01583.x/full). *Journal of Animal Ecology*, 78(6), 1113-1123.
