\name{expmFrechet}
\title{Frechet Derivative of the Matrix Exponential}
\alias{expmFrechet}
\encoding{UTF-8}
\description{
  Compute the Frechet (actually \sQuote{Fréchet}) derivative of the
  matrix exponential operator.
}
\usage{
expmFrechet(A, E, method = c("SPS", "blockEnlarge"), expm = TRUE)
}
\arguments{
  \item{A}{square matrix (\eqn{n \times n}{n x n}).}
  \item{E}{the \dQuote{small Error} matrix,
    used in \eqn{L(A,E) = f(A + E, A)}}%% FIXME
  \item{method}{string specifying the method / algorithm; the default
    \code{"SPS"} is \dQuote{Scaling + Pade + Squaring} as in the
    algorithm 6.4 below; otherwise see the \sQuote{Details} section.}
  \item{expm}{logical indicating if the matrix exponential itself, which
    is computed anyway, should be returned as well.}
}
\details{
  Calculation of \eqn{e^A} and the Exponential Frechet-Derivative
  \eqn{L(A,E)}.

  When \code{method = "SPS"} (by default), the
  with the Scaling - Padé - Squaring Method is used, in
  an R-Implementation of Al-Mohy and Higham (2009)'s Algorithm 6.4.
  \describe{
    \item{Step 1:}{Scaling (of A and E)}
    \item{Step 2:}{Padé-Approximation of \eqn{e^A} and \eqn{L(A,E)}}
    \item{Step 3:}{Squaring (reversing step 1)}
  }

  \code{method = "blockEnlarge"} uses the matrix identity of
  %% FIXME use nice LaTeX \deqn{f(\left{ .... \right} ) }
  \deqn{f([A E ; 0 A ]) = [f(A) Df(A); 0 f(A)]} for the \eqn{2n \times
    2n}{(2n) x (2n)} block matrices where \eqn{f(A) := expm(A)} and
  \eqn{Df(A) := L(A,E)}.  Note that \code{"blockEnlarge"} is much
  simpler to implement but slower (CPU time is doubled for \eqn{n = 100}).
}
\value{
  a list with components
  \item{expm}{if \code{expm} is true, the matrix exponential
    (\eqn{n \times n}{n x n} matrix).}
  \item{Lexpm}{the Exponential-Frechet-Derivative \eqn{L(A,E)}, a matrix
    of the same dimension.}
}
\references{see \code{\link{expmCond}}.}
\author{Michael Stadelmann (final polish by Martin Maechler).}
\seealso{
 \code{\link{expm.Higham08}} for the matrix exponential.
 \code{\link{expmCond}} for exponential condition number computations
  which are based on \code{expmFrechet}.
}
\examples{
(A <- cbind(1, 2:3, 5:8, c(9,1,5,3)))
E <- matrix(1e-3, 4,4)
(L.AE <- expmFrechet(A, E))
all.equal(L.AE, expmFrechet(A, E, "block"), tolerance = 1e-14) ## TRUE
}
\keyword{algebra}
\keyword{math}

