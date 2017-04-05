#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP moveHMM_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP moveHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP moveHMM_trMatrix_rcpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"moveHMM_dexp_rcpp",       (DL_FUNC) &moveHMM_dexp_rcpp,        3},
    {"moveHMM_dgamma_rcpp",     (DL_FUNC) &moveHMM_dgamma_rcpp,      3},
    {"moveHMM_dlnorm_rcpp",     (DL_FUNC) &moveHMM_dlnorm_rcpp,      3},
    {"moveHMM_dvm_rcpp",        (DL_FUNC) &moveHMM_dvm_rcpp,         3},
    {"moveHMM_dweibull_rcpp",   (DL_FUNC) &moveHMM_dweibull_rcpp,    3},
    {"moveHMM_dwrpcauchy_rcpp", (DL_FUNC) &moveHMM_dwrpcauchy_rcpp,  3},
    {"moveHMM_nLogLike_rcpp",   (DL_FUNC) &moveHMM_nLogLike_rcpp,   13},
    {"moveHMM_trMatrix_rcpp",   (DL_FUNC) &moveHMM_trMatrix_rcpp,    3},
    {NULL, NULL, 0}
};

void R_init_moveHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
