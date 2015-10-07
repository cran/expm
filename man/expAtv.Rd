\name{expAtv}
\title{Compute Matrix Exponential  exp(A t) * v  directly}
\alias{expAtv}
\alias{expmv}% - the MATLAB name "mv" = Matrix-Vector (in Stats/R mv:= multivariate)
\description{
  Compute \eqn{\exp(A t) * v} directly, without  evaluating \eqn{\exp(A)}.
}
\usage{
expAtv(A, v, t = 1,
       method = "Sidje98",
       rescaleBelow = 1e-6,
       tol = 1e-07, btol = 1e-07, m.max = 30, mxrej = 10,
       verbose = getOption("verbose"))
}
\arguments{
  \item{A}{n x n matrix}
  \item{v}{n - vector}
  \item{t}{number (scalar);}% ___ FIXME allow vector ? ___}
  \item{method}{a string indicating the method to be used; there's only
    one currently; we would like to add newer ones.}
  \item{rescaleBelow}{if \code{\link{norm}(A,"I")} is smaller than \code{rescaleBelow},
    rescale \code{A} to norm 1 and \code{t} such that \eqn{A t} remains
    unchanged.  This step is in addition to Sidje's original algorithm
    and easily seen to be necessary even in simple cases (e.g., \eqn{n = 3}).}

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

((NOT yet available!))\cr%% i.e., in expm  package, but see ../Misc/README and ../Misc/expmv/
  Al-Mohy, A. and Higham, N. (2011).
  Computing the Action of the Matrix Exponential, with an Application
  to Exponential Integrators.
  \emph{SIAM Journal on Scientific Computing}, \bold{33}(2), 488--511.
%% @article{doi:10.1137/100788860,
%%   author =	{Al-Mohy, A. and Higham, N.},
%%   title =	{Computing the Action of the Matrix Exponential, with an
%%		 Application to Exponential Integrators},
%%   journal =	{SIAM Journal on Scientific Computing},
%%   volume =	33,
%%   number =	2,
%%   pages =	{488-511},
%%   year =	2011,
%%   doi =		{10.1137/100788860},
%%   URL =		{http://epubs.siam.org/doi/abs/10.1137/100788860},
%% }
}
\author{Ravi Varadhan, Johns Hopkins University;
  Martin Maechler (cosmetic, generalization to sparse matrices;
  rescaling (see \code{rescaleBelow}).
}
% \note{
% }

\seealso{
  \code{\link{expm}}%,....
}
\examples{%% tests are at end of ../tests/ex2.R
source(system.file("demo", "exact-fn.R", package = "expm"))
##-> rnilMat() ; xct10
set.seed(1)
(s5 <- Matrix(m5 <- rnilMat(5))); v <- c(1,6:9)
(em5 <- expm(m5))
r5  <- expAtv(m5, v)
r5. <- expAtv(s5, v)
stopifnot(all.equal(r5, r5., tolerance = 1e-14),
          all.equal(c(em5 \%*\% v), r5$eAtv))

v <- 10:1
with(xct10, all.equal(expm(m), expm))
all.equal(c(xct10$expm \%*\% v),
          expAtv(xct10$m, v)$eAtv)
}
\keyword{algebra}
\keyword{math}
