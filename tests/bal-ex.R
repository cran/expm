library(expm)
source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)## -> assertError()...

demo(balanceTst) #-> the function definition and the first few examples

dm4. <- dgebal(m4)
storage.mode(m4) <- "integer"
stopifnot(identical(dm4., dgebal(m4)))

expm(m)
expm(m,"Pade") ## are different indeed {when bug still existed}
expm(m,"R_Pade")# same as Pade


## a non-empty ``non-balanced'' example  ---

expm.t.identity(m4, "Ward")

m6 <- zeroTrace(matrix(outer(2^(-8:9),c(-1,1)), 6,6)); m6
m6[lower.tri(m6)] <- 0 ## plus one non-zero
m6[4,2] <- 77
p <- c(6,4,5,2:1,3); m6 <- m6[p,p]
expm.t.identity(m6, "Ward") ##  difference; indeed
expm(m6) # is very different from
expm(m6,"R_Pade")

str(dm6 <- dgebalTst(m6))
## Now, that's interesting:
##
## 1.  'S' scales *more* (2 .. 5) than just (2:4 == i1:i2) !
##
## 2.  'B' has quite different scaling and it does (must!) obey rule
##     scale i1:i2 only
##
## 3. 'B'(oth) is better than  "P" and "S" separately:
##
kappa(eigen(m6)$vectors)#      597.5588
kappa(eigen(dm6$P$z)$vectors)# 597.5588
kappa(eigen(dm6$S$z)$vectors)#  42.58396
kappa(eigen(dm6$B$z)$vectors)#  22.20266


## An n=17 example where octave's expm() is wrong too
m17 <- matrix(c(10,0, 0, 2, 3,-1, 0, 0, 0, 0, 0, 4, 0, 5, 0, 0,-2,
                0, 0, 0, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 7, 0,
                0, 0,10, 0, 0,-4, 9, 0, 0, 0,-5, 0,-6, 0, 0, 0, 0,
                0, 0,-7, 0, 0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0,11, 0,
                0, 0, 0, 0, 0, 0,12, 0, 0, 0, 0, 0,-8, 0, 0, 0, 0,
                0, 0,-9, 0, 0, 0, 0, 0, 0,-10,0,13,14,-11,-12,-13, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0,10, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0,-14,16,0,-10,0,17, 0, 0, 0, 0, 0, 0,
                0, 0,-16,0, 0,18,19, 0, 0, 0, 0, 0, 0, 0,20, 0, 21,
                22,0, 0, 0, 0, 0,-17,0, 0, 0,-10,-19,-20,0,0,0, 0,
                0,-21,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0,23,24, 0,25,26, 0, 0,27,-22,0,28,-23,0,-24,
                0,-25,0,29, 0, 0, 0, 0, 0, 0, 0,30,31, 0, 0, 0, 0,
                0, 0,-26,32,0, 0, 0, 0, 0,-27,0,33,34, 0, 0, 0, 0,
                0,-28,-29,0,0, 0,35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0,36,37, 0, 0, 0, 0, 0, 0, 0, 0,-10),
              17, 17)
str(dm17 <- dgebalTst(m17))
sapply(dm17[1:3], `[[`, "scale")

## The balancing was really rather harmful -- cond(V) *not* improved:
condX <- function(x) kappa(x, exact=TRUE)
condX(eigen(m17)$vectors)#      8.9e16
condX(eigen(dm17$P$z)$vectors)# 1.37e17
condX(eigen(dm17$S$z)$vectors)# 1.44e17
condX(eigen(dm17$B$z)$vectors)# 1.43e17 (very slightly smaller)
