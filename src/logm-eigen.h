/*  ===== File part of R package expm =====
 *
 *  logm-eigen.h
 *
 *  Created by Christophe Dutang on 13/05/08.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>



#include "R_NLS_locale.h"
//#include "logm.h"

SEXP do_logm_eigen(SEXP x, SEXP tolin);
void logm_eigen(double *x, int n, double *z, double tol);


