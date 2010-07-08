
#ifndef R_PKG_EXPM_H
#define R_PKG_EXPM_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "locale.h"

typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;

SEXP do_expm(SEXP x, SEXP kind);
void expm(double *x, int n, double *z, precond_type precond_kind);

#endif /* R_PKG_EXPM_H */
