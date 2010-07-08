#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

#include "locale.h"

/* The C API :*/
void matpow(double *x, int n, int k, double *z);

/* as .Call()ed from R */
SEXP R_matpow(SEXP x, SEXP k);

