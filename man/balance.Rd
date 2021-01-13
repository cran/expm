\name{balance}
\alias{balance}
\alias{dgebal}
\title{Balance a Square Matrix via LAPACK's DGEBAL}
\description{
  Balance a square matrix via LAPACK's \code{DGEBAL}.
  This is an \R interface, mainly used for experimentation.

  This LAPACK routine is used internally for Eigenvalue decompositions,
  but also, in Ward(1977)'s algorithm for the matrix exponential.

  The name \code{balance()} is preferred nowadays, and \dQuote{dgebal()}
  has been deprecated (finally, after 9 years ...).
}
\usage{
balance(A, job = c("B", "N", "P", "S"))
## Deprecated now:
##   dgebal(A, job = c("B", "N", "P", "S"))
}
\arguments{
  \item{A}{a square (\eqn{n\times n}{n x n}) numeric matrix.}
  \item{job}{a one-letter string specifying the \sQuote{job} for DGEBAL.
    \describe{
      \item{P}{\bold{P}ermutation}
      \item{S}{\bold{S}caling}
      \item{B}{\bold{B}oth permutation and scaling}
      \item{N}{\bold{N}one}
    }
  }
}
\value{
  A list with components
  \item{z}{the transformation of matrix \code{A}, after permutation and
    or scaling.}
  \item{scale}{numeric vector of length \eqn{n}, containing the
    permutation and/or scale vectors applied.}
  \item{i1,i2}{integers (length 1) in \eqn{\{1,2,\dots,n\}}, denoted by
    \code{ILO} and \code{IHI} respectively in the LAPACK
    documentation.  Only relevant for \code{"P"} or \code{"B"}, they
    describe where permutations and where scaling took place; see the
    Details section.
  }
}
\details{
  An excerpt of the LAPACK documentation about DGEBAL(), describing the
  result
  \describe{
    \item{i1 ("ILO")}{(output) integer}
    \item{i2 ("IHI")}{(output) integer\cr
      \code{i1} and \code{i2} are set to integers such that on exit
      \code{z[i,j] = 0} if \eqn{i > j} and \eqn{j = 1,...,i1-1} or \eqn{i = i2+1,...,n}.

      If \code{job = 'N'} or \code{'S'}, \code{i1 = 1} and \code{i2 = n}.}

    \item{scale}{(output) numeric vector of length \code{n}.
      Details of the permutations and scaling factors applied to
      \code{A}.  If \code{P[j]} is the index of the row and column interchanged
      with row and column \code{j} and \code{D[j]} is the scaling factor
      applied to row and column j, then
      \code{scale[j] = P[j]} for \eqn{j = 1,...,i1-1}\cr
      \code{         = D[j]} for \eqn{j = i1,...,i2},\cr
      \code{         = P[j]} for \eqn{j = i2+1,...,n}.

      The order in which the interchanges are made is \code{n} to \code{i2+1},
      then \code{1} to \code{i1-1}.}
  }
  Look at the LAPACK documentation for more details.
}
\references{
  LAPACK Reference Manual
}
\author{Martin Maechler}
\seealso{\code{\link{eigen}}, \code{\link{expm}}.
}
\examples{
m4 <- rbind(c(-1,-1, 0, 0),
            c( 0, 0,10,10),
            c( 0, 0,10, 0),
            c( 0,10, 0, 0))
(b4 <- balance(m4))

## --- for testing and didactical reasons : ----

demo(balanceTst) # also defines the balanceTst() function
                # which in its tests ``defines'' what
                # the return value means, notably (i1,i2,scale)
}
\keyword{array}
\keyword{arith}
