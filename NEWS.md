
# moveHMM 1.12

* Change `arma::is_finite()` to `std::isfinite()` after Armadillo update
* Require >= 3 observations per track in `prepData` (to have one turning angle)

# moveHMM 1.11

* Add built-in functions for von Mises and wrapped Cauchy distributions (to remove dependence on CircStats, as requested by CRAN): `dvm`, `rvm`, `dwrpcauchy`, `rwrpcauchy`.

# moveHMM 1.10 

* Add `splitAtGaps` function
