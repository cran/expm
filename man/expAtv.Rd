\name{expAtv}
\alias{expAtv}
\title{Compute Matrix Exponential  exp(A t) * v  directly}
\description{Performs exp(A t) \%*\% v  directly  w/o  evaluating exp(A)
%%  Originally by

}
\usage{
expAtv(A, v, t = 1,
       method = "Sidje98",
       tol = 1e-07, btol = 1e-07, m.max = 30, mxrej = 10,
       verbose = getOption("verbose"))
}
\arguments{
  \item{A}{n x n matrix}
  \item{v}{n - vector}
  \item{t}{number (scalar);}% ___ FIXME allow vector ? ___}
  \item{method}{a string indicating the method to be used; there's only
    one currently.}
  \item{tol, btol}{tolerances; these are tuning constants of the
    "Sidje1998" method which the user should typically \emph{not} change.}
  \item{m.max, mxrej}{integer constants you should only change if you know what you're doing}
  \item{verbose}{flag indicating if the algorithm should be verbose..}
}
\value{a list with components
  \item{eAtv}{.....fixme...}
}
% \details{

% }
\references{
  Roger B. Sidje (1998)
  EXPOKIT: Software Package for Computing Matrix Exponentials.
  \emph{ACM - Transactions On Mathematical Software} \bold{24}(1), 130--156.
}
\author{Ravi Varadhan, Johns Hopkins University;
  Martin Maechler (cosmetic, generalization to sparse matrices).
}
% \note{
% }

\seealso{
  \code{\link{expm}},....
}
\examples{
source(system.file("demo", "exact-fn.R", package = "expm"))
##-> rnilMat() ; xct10
set.seed(1)
(s5 <- Matrix(m5 <- rnilMat(5))); v <- c(1,6:9)
(em5 <- expm(m5))
r5  <- expAtv(m5, v)
r5. <- expAtv(s5, v)
stopifnot(all.equal(r5, r5., tol = 1e-14),
          all.equal(c(em5 \%*\% v), r5$eAtv))

v <- 10:1
with(xct10, all.equal(expm(m), expm))
all.equal(c(xct10$expm \%*\% v),
          expAtv(xct10$m, v)$eAtv)

(doExtras <- interactive() || nzchar(Sys.getenv("R_expm_check_extra")) ||
 identical("true", unname(Sys.getenv("R_MM_PKG_CHECKING"))))

n <- 500
A <- bandSparse(n,n, -1:1, diag = list(-(2:n), -5*(1:n), 1:(n-1)))
v <- 100*(500:1)
t.v <- system.time(rr  <- expAtv(A, v=v))
if(doExtras) { ## this is an order of magnitude slower :
 t.A <- system.time(w. <- (eA <- expm(A, "Higham08")) \%*\% v)
 stopifnot(all.equal(rr$eAtv, as.numeric(w.)))
 print( mean((t.A / t.v)[c(1,3)]) )## 23.57  {nb-mm3}
}% if(doExtras)
}
\keyword{algebra}
\keyword{math}
