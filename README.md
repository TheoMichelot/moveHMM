# moveHMM
An R package for animal movement modelling using hidden Markov models.

## Installation instructions

1. Install package dependencies : CircStats, sp, Rcpp, RcppArmadillo (from CRAN, using `install.packages` for example).

2. Install the package devtools from CRAN, and load it.

3. Install moveHMM from the GitHub repository using `install.github("TheoMichelot/moveHMM")`.

Alternatively, on Linux, you can :

Download ZIP archive of the repository (moveHMM-master.zip). Extract this archive to moveHMM-master. 
Then, build the package from the command line using `R CMD build moveHMM-master`. This will create a tarbal of the package. 
Finally, open R and install the package using `install.packages`, e.g. `install.packages("moveHMM_0.1.tar.gz")`.
