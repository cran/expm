\name{logm}
\alias{logm}
\title{Matrix Logarithm}
\description{
  This function computes the (principal) matrix logarithm of a square matrix.
  A logarithm of a matrix \eqn{A} is \eqn{L} such that \eqn{A= e^L}
  (meaning \code{A == expm(L)}), see the documentation for the matrix
  exponential, \code{\link{expm}}, which can be defined
  as \deqn{e^L := \sum_{r=0}^\infty L^r/r! .}{sum[r = 0,..,Inf; L^r/r!] .}
}
\usage{
logm(x, method = c("Higham08", "Eigen"),
%      order = 8, trySym = TRUE,
     tol = .Machine$double.eps)
}
\arguments{
  \item{x}{a square matrix.}
  \item{method}{a string specifying the algorithmic method to be used.
    The default uses the algorithm by Higham(2008).

    The simple \code{"Eigen"} method tries to diagonalise the matrix
    \code{x}; if that is not possible, it raises an error.
  }
%   \item{order}{only for \code{method = "Eigen"} .. FIXME}
%   \item{trySym}{only for \code{method = "Eigen"} .. FIXME}
  \item{tol}{a given tolerance used to check if \code{x} is
    computationally singular when \code{method = "Eigen"}.}
}
\details{
  The exponential of a matrix is defined as the infinite Taylor series
  \deqn{e^M = \sum_{k = 1}^\infty \frac{M^k}{k!}.}{%
    exp(M) = I + M + M^2/2! + M^3/3! + \dots.}
  The matrix logarithm of \eqn{A} is a matrix \eqn{M} such that
  \eqn{exp(M) = A}. Note that there typically are an infinite number
  number of such matrices, and we compute the \emph{prinicipal} matrix
  logarithm, see the references.

  Method \code{"Higham08"} works via \dQuote{\emph{inverse} scaling and
    squaring}, and from the Schur decomposition, applying a matrix
  square root computation.  It is somewhat slow but also works for
  non-diagonalizable matrices.
}
\references{
  Higham, N.~J. (2008).
  \emph{Functions of Matrices: Theory and Computation};
  Society for Industrial and Applied Mathematics, Philadelphia, PA, USA.

  The Matrix Logarithm is very nicely defined by Wikipedia,
  \url{http://en.wikipedia.org/wiki/Matrix_logarithm}.
}
\value{
  A matrix \sQuote{as \code{x}} with the matrix logarithm of \code{x},
  i.e., \code{all.equal( expm(logm(x)), x, tol)} is typically true for
  quite small tolerance \code{tol}.
}
\seealso{
  \code{\link{expm}}
}
\author{
  Method \code{"Higham08"} was implemented by Michael Stadelmann as part of his
  master thesis in mathematics, at ETH Zurich;
  the \code{"Eigen"} method by Christophe Dutang.
}
\examples{
m <- diag(2)
logm(m)
expm(logm(m))

## Here, logm() is barely defined, and Higham08 has needed an amendment
## in order for not to loop forever:
D0 <- diag(x=c(1, 0.))
(L. <- logm(D0))
stopifnot( all.equal(D0, expm(L.)) )

## A matrix for which clearly no logm(.) exists:
(m <- cbind(1:2, 1))
(l.m <- try(logm(m))) ## all NA {Warning in sqrt(S[ij, ij]) : NaNs produced}
## on r-patched-solaris-x86, additionally gives
##    Error in solve.default(X[ii, ii] + X[ij, ij], S[ii, ij] - sumU) :
##     system is computationally singular: reciprocal condition number = 0
##    Calls: logm ... logm.Higham08 -> rootS -> solve -> solve -> solve.default
if(!inherits(l.m, "try-error")) stopifnot(is.na(l.m))
## The "Eigen" method  ``works''  but wrongly :
expm(logm(m, "Eigen"))
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
