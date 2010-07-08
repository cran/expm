\name{expm}
\alias{expm}
\alias{mexp}% for now
\title{Matrix Exponential}
\description{
  This function computes the exponential of a square matrix
  \eqn{A}, defined as the sum from \eqn{r=0} to infinity of
  \eqn{A^r/r!}.
  Several methods are provided. The Taylor series and Padé
  approximation are very importantly combined with scaling and squaring.
}
\usage{
expm(x, method = c("Higham08.b", "Higham08",
                   "Ward77", "Pade", "Taylor", "PadeO", "TaylorO",
                   "R_Eigen", "R_Pade", "R_Ward77", "hybrid_Eigen_Ward"),
     order = 8, trySym = TRUE, tol = .Machine$double.eps,
     preconditioning = c("2bal", "1bal", "buggy"))
}
\arguments{
  \item{x}{a square matrix.}
  \item{method}{\code{"Higham08.b"}, \code{"Ward77"}, \code{"Pade"} or
    \code{"Taylor"}, etc;
    The default is now \code{"Higham08.b"} which uses Higham's 2008
    algorithm with \emph{additional} balancing preconditioning,
    see \code{\link{expm.Higham08}}.

    The versions with "*O" call the
    original Fortran code, whereas the first ones call the BLAS-using
    and partly simplified newer code.\cr
    \code{"R_Pade"} uses an \R-code version of \code{"Pade"} for
    didactical reasons, and\cr
    \code{"R_Ward77"} uses an \R version of \code{"Ward77"}, still based
    on LAPACK's \code{dgebal}, see \R interface \code{\link{dgebal}}.
    This has enabled us to diagnose and fix the bug in the original
    octave implementation of \code{"Ward77"}.
    \code{"R_Eigen"} tries to diagonalise the matrix \code{x}, if not possible,
    \code{"R_Eigen"} raises an error. \code{"hybrid_Eigen_Ward"} method also
    tries to diagonalise the matrix \code{x}, if not possible, it uses
    \code{"Ward77"} algorithm.
    }
  \item{order}{an integer, the order of approximation to be used, for
    the "Pade" and "Taylor" methods.  The best value for this depends on
    machine precision (and slightly on \code{x}) but for the current
    double precision arithmetic, one recommendation %% and Matlab !!!
    uses \code{order = 6} unconditionally;\cr
    %% NOTA BENE:  Matlab uses order = 6  !!!
    our default, \code{8}, is from Ward(1977, p.606)'s recommendation.
  }
  \item{trySym}{logical indicating if \code{method = "R_Eigen"} should use
    \code{\link{isSymmetric}(x)} and take advantage for (almost)
    symmetric matrices.}
  \item{tol}{a given tolerance used to check if \code{x} is
    computationally singular when \code{method = "hybrid_Eigen_Ward"}.}
  \item{preconditioning}{a string specifying which implementation of
    Ward(1977) should be used when \code{method = "Ward77"}.
    %%% FIXME explain more
  }
}
\details{
  The exponential of a matrix is defined as the infinite Taylor series
  \deqn{e^M = \sum_{k = 1}^\infty \frac{M^k}{k!}.}{%
    exp(M) = I + M + M^2/2! + M^3/3! + \dots}

  For the "Pade" and "Taylor" methods, there is an \code{"accuracy"}
  attribute of the result. It is an upper bound for the L2 norm of the
  Cauchy error \code{expm(x, *, order + 10) - expm(x, *, order)}.
}
\value{
  The matrix exponential of \code{x}.
}
\seealso{
  The package vignette for details on the algorithms and calling the
  function from external packages.
}
\author{
  The \code{"Ward77"} method:\cr
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}, and Christophe
  Dutang, based on code translated by Doug Bates and Martin Maechler
  from the implementation of the corresponding Octave function
  contributed by A. Scottedward Hodel \email{A.S.Hodel@eng.auburn.edu}.

  The \code{"PadeO"} and  \code{"TaylorO"} methods:\cr
  Marina Shapira (U Oxford, UK) and David Firth (U Warwick, UK);

  The \code{"Pade"} and  \code{"Taylor"} methods are slight
  modifications of the "*O" ([O]riginal versions) methods,
  by Martin Maechler, using BLAS and LINPACK where possible.

  The \code{"hybrid_Eigen_Ward"} method by Christophe Dutang is a C
  translation of \code{"R_Eigen"} method by Martin Maechler.

  The \code{"Higham08"} and \code{"Higham08.b"} (current default) were
  written by Michael Stadelmann, see \code{\link{expm.Higham08}}.
}
\references{
  Ward, R. C. (1977).  Numerical computation
  of the matrix exponential with accuracy estimate.
  \emph{SIAM J. Num. Anal.} \bold{14}, 600--610.

  Moler, C and van Loan, C (2003).  Nineteen dubious ways to compute
  the exponential of a matrix, twenty-five years later.
  \emph{SIAM Review} \bold{45}, 3--49.  At
  \url{http://epubs.siam.org/sam-bin/dbq/article/41801}
}
\note{
  For a good general discussion of the matrix exponential
  problem, see Moler and van Loan (2003).
}
\examples{
x <- matrix(c(-49, -64, 24, 31), 2, 2)
expm(x)

## ----------------------------
## Test case 1 from Ward (1977)
## ----------------------------
test1 <- t(matrix(c(
    4, 2, 0,
    1, 4, 1,
    1, 1, 4), 3, 3))
expm(test1, method="Pade")
## Results on Power Mac G3 under Mac OS 10.2.8
##                    [,1]               [,2]               [,3]
## [1,] 147.86662244637000 183.76513864636857  71.79703239999643
## [2,] 127.78108552318250 183.76513864636877  91.88256932318409
## [3,] 127.78108552318204 163.67960172318047 111.96810624637124
## -- these agree with ward (1977, p608)

## Compare with the naive "R_Eigen" method:
try(
expm(test1, method="R_Eigen")
) ## now gives an error from solve !
##
## older result was
##                   [,1]                [,2]               [,3]
##[1,] 147.86662244637003  88.500223574029647 103.39983337000028
##[2,] 127.78108552318220 117.345806155250600  90.70416537273444
##[3,] 127.78108552318226  90.384173332156763 117.66579819582827
## -- hopelessly inaccurate in all but the first column.
##
##
## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
test2 <- t(matrix(c(
    29.87942128909879, .7815750847907159, -2.289519314033932,
    .7815750847907159, 25.72656945571064,  8.680737820540137,
   -2.289519314033932, 8.680737820540137,  34.39400925519054),
           3, 3))
expm(test2, method="Pade")
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692357 -18231880972009844 -30475770808580828
##[2,] -18231880972009852  60605228702227024 101291842930256144
##[3,] -30475770808580840 101291842930256144 169294411240859072
## -- which agrees with Ward (1977) to 13 significant figures
expm(test2, method="R_Eigen")
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692405 -18231880972009100 -30475770808580196
##[2,] -18231880972009160  60605228702221760 101291842930249376
##[3,] -30475770808580244 101291842930249200 169294411240850880
## -- in this case a very similar degree of accuracy.
##
## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
test3 <- t(matrix(c(
    -131, 19, 18,
    -390, 56, 54,
    -387, 57, 52), 3, 3))
expm(test3, method="Pade")
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608.
expm(test3, method="R_Eigen")
##                   [,1]               [,2]                [,3]
##[1,] -1.509644158796182 0.3678794391103086 0.13533528117547022
##[2,] -5.632570799902948 1.4715177585023838 0.40600584352641989
##[3,] -4.934938326098410 1.1036383173309319 0.54134112676302582
## -- in this case, a similar level of agreement with Ward (1977).
##
## ----------------------------
## Test case 4 from Ward (1977)
## ----------------------------
test4 <-
    structure(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10,
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    .Dim = c(10, 10))
attributes(expm(test4, method="Pade"))
max(abs(expm(test4, method="Pade") - expm(test4, method="R_Eigen")))
##[1] 8.746826694186494e-08
## -- here mexp2 is accurate only to 7 d.p., whereas mexp
##    is correct to at least 14 d.p.
##
## Note that these results are achieved with the default
## settings order=8, method="Pade" -- accuracy could
## presumably be improved still further by some tuning
## of these settings.

##
## example of computationally singular matrix
##
m <- matrix(c(0,1,0,0), 2,2)
try(
expm(m, m="R_Eigen")
)
## error since m is computationally singular
expm(m, m="hybrid")
## hybrid use the Ward77 method
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
