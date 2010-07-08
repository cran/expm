dgebal <- function(A, job = c("B","P","S"))
    .Call("R_dgebal", A, match.arg(job))

