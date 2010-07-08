\name{sqrtm}
\alias{sqrtm}
\title{Matrix Square Root}
\description{
  This function computes the matrix square root of a square matrix.
  The sqrt of a matrix \eqn{A} is \eqn{S} such that \eqn{A = S S}.
}
\usage{
sqrtm(x)
}
\arguments{
  \item{x}{a square matrix.}
}
\details{
  The matrix square root \eqn{S} of \eqn{M}, \eqn{S = sqrtm(M)} is
  defined as one (the \dQuote{principal}) \eqn{S} such that
  \eqn{S S = S^2 = M}, (in \R, \code{all.equal( S \%*\% S , M )}).

  The method works from the Schur decomposition.
}
\value{
  A matrix \sQuote{as \code{x}} with the matrix sqrt of \code{x}.
}
\references{
  Higham, N.~J. (2008).
  \emph{Functions of Matrices: Theory and Computation};
  Society for Industrial and Applied Mathematics, Philadelphia, PA, USA.
}
\seealso{
  \code{\link{expm}}, \code{\link{logm}}
}
\author{
  Michael Stadelmann wrote the first version.
}
\examples{
m <- diag(2)
sqrtm(m) == m # TRUE

(m <- rbind(cbind(1, diag(1:3)),2))
sm <- sqrtm(m)
sm
zapsmall(sm \%*\% sm) # Zap entries ~= 2e-16
stopifnot(all.equal(m, sm \%*\% sm))
}
\keyword{algebra}
\keyword{math}
\encoding{UTF-8}
