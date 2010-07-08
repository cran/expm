\name{expmCond}
\title{Exponential Condition Number of a Matrix}
\Rdversion{1.1}
\alias{expmCond}
\description{
  Compute the exponential condition number of a matrix, either with
  approximation methods, or exactly and very slowly.
}
\usage{
expmCond(A, method = c("1.est", "F.est", "exact"),
         expm = TRUE, abstol = 0.1, reltol = 1e-6,
         give.exact = c("both", "1.norm", "F.norm"))
}
\arguments{
  \item{A}{a square matrix}
  \item{method}{a string; either compute 1-norm or F-norm
    \emph{approximations}, or compte these \emph{exactly}.}
  \item{expm}{logical indicating if the matrix exponential itself, which
    is computed anyway, should be returned as well.}
  \item{abstol, reltol}{for \code{method = "F.est"}, numerical \eqn{\ge 0}{>= 0},
    as \emph{absolute} and \emph{relative} error tolerance.}
  \item{give.exact}{for \code{method = "exact"}, specify if only the
    1-norm, the Frobenius norm, or both are to be computed.}
}
\details{
  \code{method = "exact"}, aka Kronecker-Sylvester algorithm, computes a
  Kronecker matrix of dimension \eqn{n^2 \times n^2}{n^2 x n^2} and
  hence, with \eqn{O(n^5)} complexity, is prohibitely slow for
  non-small \eqn{n}.  It computes the \emph{exact} exponential-condition
  numbers for both the Frobenius and/or the 1-norm, depending on
  \code{give.exact}.

  The two other methods compute approximations, to these norms, i.e.,
  \bold{est}imate them, using algorithms from Higham, chapt.~3.4, both
  with complexity \eqn{O(n^3)}.

  %% FIXME:  Say more
}
\value{
  when \code{expm = TRUE}, for \code{method = "exact"}, a
  \code{\link{list}} with components
  \item{expm}{containing the matrix exponential, \code{\link{expm.Higham08}(A)}.}
  \item{expmCond(F|1)}{numeric scalar, (an approximation to) the (matrix
    exponential) condition number, for either the 1-norm
    (\code{expmCond1}) or the Frobenius-norm (\code{expmCondF}).}

  When \code{expm} is false and \code{method} one of the approximations
  (\code{"*.est"}), the condition number is returned directly (i.e.,
  \code{\link{numeric}} of length one).
}
\references{
  Awad H. Al-Mohy and Nicholas J. Higham (2009).
  \emph{Computing Fréchet Derivative of the Matrix Exponential, with an application
    to Condition Number Estimation}; MIMS EPrint 2008.26; Manchester
  Institute for Mathematical Sciences, U. Manchester, UK.
  \url{http://eprints.ma.man.ac.uk/1218/01/covered/MIMS_ep2008_26.pdf}

  Higham, N.~J. (2008).
  \emph{Functions of Matrices: Theory and Computation};
  Society for Industrial and Applied Mathematics, Philadelphia, PA, USA.

  Michael Stadelmann (2009) \emph{Matrixfunktionen} ...
  Master's thesis; see reference in \code{\link{expm.Higham08}}.
}
\author{
  Michael Stadelmann (final polish by Martin Maechler).
}
\seealso{
 \code{\link{expm.Higham08}} for the matrix exponential.
}
\examples{
set.seed(101)
(A <- matrix(round(rnorm(3^2),1), 3,3))

eA <- expm.Higham08(A)
stopifnot(all.equal(eA, expm::expm(A), tol= 1e-15))

C1 <- expmCond(A, "exact")
C2 <- expmCond(A, "1.est")
C3 <- expmCond(A, "F.est")
all.equal(C1$expmCond1, C2$expmCond, tol= 1e-15)# TRUE
all.equal(C1$expmCondF, C3$expmCond)# relative difference of 0.001...
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
