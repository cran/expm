\name{matStig}
\alias{matStig}
\docType{data}
\title{Stig's "infamous" Example Matrix}
\description{
  Stig Mortensen wrote on Oct 22, 2007 to the authors of the \pkg{Matrix}
  package with subject \dQuote{Strange result from expm}.
  There, he presented the following \eqn{8 \times 8}{8 * 8} matrix for
  which the Matrix \code{expm()} gave a \dQuote{strange} result.
  As we later researched, the result indeed was wrong: the correct
  entries were wrongly permuted.  The reason has been in the underlying
  source code in Octave from which it had been ported to \pkg{Matrix}.
}
\usage{data(matStig)}
\author{Martin Maechler}
\examples{
data(matStig)

as(matStig, "sparseMatrix") # since that prints more nicely.

## For more compact printing:
op <- options(digits = 4)

E1 <- expm(matStig, "Ward77", preconditioning="buggy") # the wrong result
as(E1, "sparseMatrix")
str(E2 <- expm(matStig, "Pade"))# the correct one (has "accuracy" attribute)
as(E2, "sparseMatrix")
attr(E2,"accuracy") <- NULL   # don't want it below
E3 <- expm(matStig, "R_Eigen")  # even that is fine here
all.equal(E1,E2) # not at all equal (rel.difference >~= 1.)
stopifnot(all.equal(E3,E2)) # ==

##________ The "proof" that "Ward77" is wrong _________
M <- matStig
Et1 <- expm(t(M), "Ward77", precond= "buggy")
Et2 <- expm(t(M), "Pade"); attr(Et2,"accuracy") <- NULL
all.equal(Et1, t(E1)) # completely different (rel.diff ~ 1.7 (platform dep.))
stopifnot(all.equal(Et2, t(E2))) # the same (up to tolerance)

options(op)
}
\keyword{array}
\keyword{datasets}
