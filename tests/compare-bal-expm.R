library(expm)

source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)

## Compare the two different "balance" pre-conditioning versions in Ward77:

set.seed(1)
mList <- lapply(integer(100), function(...) rSpMatrix(20, nnz=80))
re20 <- sapply(mList, function(M)
               relErr(expm(M, precond = "2bal"),
                      expm(M, precond = "1bal")))
re20 ## ahh.: zero or ~ 1e-13 ... good
table(re20 == 0)
summary(re20[re20 != 0])
## Pentium M (ubuntu)
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
## 2.593e-14 8.703e-14 1.282e-13 2.434e-13 4.177e-13 6.295e-13
