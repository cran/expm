/*  ===== File part of R package expm =====
 *
 *  expm-eigen.h
 *
 *  Created by Christophe Dutang on 27/02/08.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>



#include "R_NLS_locale.h"
#include "expm.h"

SEXP do_expm_eigen(SEXP x, SEXP tolin);
void expm_eigen(double *x, int n, double *z, double tol);


