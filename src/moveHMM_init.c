#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton

/* .Call calls */
extern SEXP _moveHMM_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP _moveHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _moveHMM_trMatrix_rcpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_moveHMM_dexp_rcpp",       (DL_FUNC) &_moveHMM_dexp_rcpp,        3},
    {"_moveHMM_dgamma_rcpp",     (DL_FUNC) &_moveHMM_dgamma_rcpp,      3},
    {"_moveHMM_dlnorm_rcpp",     (DL_FUNC) &_moveHMM_dlnorm_rcpp,      3},
    {"_moveHMM_dvm_rcpp",        (DL_FUNC) &_moveHMM_dvm_rcpp,         3},
    {"_moveHMM_dweibull_rcpp",   (DL_FUNC) &_moveHMM_dweibull_rcpp,    3},
    {"_moveHMM_dwrpcauchy_rcpp", (DL_FUNC) &_moveHMM_dwrpcauchy_rcpp,  3},
    {"_moveHMM_nLogLike_rcpp",   (DL_FUNC) &_moveHMM_nLogLike_rcpp,   13},
    {"_moveHMM_trMatrix_rcpp",   (DL_FUNC) &_moveHMM_trMatrix_rcpp,    3},
    {NULL, NULL, 0}
};

void R_init_moveHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
