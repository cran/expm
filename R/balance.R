## NOTA BENE: In Matlab, there's the function balance(.) which
## calls LAPACK's  dgebal  *AND* which has the option to also return the
## transformation *diagonal* matrix D , not just the transformed matrix.
balance <- function(A, job = c("B","N", "P","S")) {
    .Call(R_dgebal,
          if(is.numeric(A)) A else as(A, "matrix"),
          match.arg(job))
}

## dgebal <- balance ## till
## 2020-07-21: Finally deprecated:
dgebal <- function(A, job = c("B","N", "P","S")) {
    .Deprecated("balance")
    .Call("R_dgebal", A, match.arg(job))
}

## Not exported, used to make  'R CMD check <pkg>'  be faster *or* more extensive:
doExtras <- function(int = interactive()) {
    int || nzchar(Sys.getenv("R_expm_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}
