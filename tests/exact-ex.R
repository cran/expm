#### Examples where we know the result "exactly"

library(expm)

source(system.file("test-tools.R", package = "expm"))## -> assertError(), rMat()
source(system.file("demo", "exact-fn.R", package = "expm"))

set.seed(321)
re <- replicate(1000,
                c(re.nilA3(rlnorm(3),function(x)expm(x,"Pade")),
                  re.nilA3(rnorm(3), function(x)expm(x,"Pade"))))

summary(t(re))
stopifnot(rowMeans(re) < 1e-15,
          apply(re, 1, quantile, 0.80) < 1e-16,
          apply(re, 1, quantile, 0.90) < 2e-15,
          apply(re, 1, max) < c(4e-14, 6e-15))

cat('Time elapsed: ', (p1 <- proc.time()),'\n') # for ``statistical reasons''


## Check *many* random nilpotent matrices:
set.seed(321)
RE <- replicate(1000,
                c(re.nilA3(rlnorm(3), function(x) expm(x, "Ward77")),
                  re.nilA3(rnorm(3),  function(x) expm(x, "Ward77"))))
stopifnot(rowMeans(RE) < 1e-15,
          apply(RE, 1, quantile, 0.80) < 1e-16,
          apply(RE, 1, quantile, 0.90) < 2e-15,
          apply(RE, 1, max) < c(4e-14, 6e-15))

print(summary(t(RE)))
epsC <- .Machine$double.eps
cat("relErr(expm(.,Pade)) - relErr(expm(.,Ward77))  in Machine_eps units:\n")
print(summary(c(re - RE)) / epsC)
##       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
## -0.6183442  0.0000000  0.0000000  1.3650410  0.1399719 94.9809161

cat('Time elapsed: ',(p2 <- proc.time())-p1,'\n') # for ``statistical reasons''

###--- A second group --- where we know the diagonalization of A ---

if(!require("Matrix"))
    q('no')
##  ------  the rest really uses 'Matrix'
##---> now use  expm::expm()  since Matrix has its own may mask the expm one
##              ^^^^^^^^^^^^

## rMat() relies on Matrix::rcond():
## Now with the change default rcondMin, this "works"
system.time(R40 <- rMat(40))
system.time(R80 <- rMat(80))

expm.safe.Eigen <- function(x, silent = FALSE) {
    r <- try(expm::expm(x, "R_Eigen"), silent = silent)
    if(inherits(r, "try-error")) NA else r
}

expmList <-
    list(Ward  = function(x) expm::expm(x, "Ward77"),
	 s.P.s = function(x) expm::expm(x, "Pade"),
	 s.P.sO= function(x) expm::expm(x, "PadeO"),
	 sPs.H08.= function(x) expm::expm.Higham08(x, balancing=FALSE),
	 sPs.H08b= function(x) expm::expm.Higham08(x, balancing= TRUE),
	 s.T.s = function(x) expm::expm(x, "Taylor"),
	 s.T.sO= function(x) expm::expm(x, "TaylorO"),
	 Eigen = expm.safe.Eigen,
	 hybrid= function(x) expm::expm(x, "hybrid")
	 )
expmL.wo.E <- expmList[names(expmList) != "R_Eigen"]


set.seed(12)
re.facMat(20, expmList)
fRE <- replicate(100, re.facMat(20, expmList))

## Now look at that:
boxplot(t(fRE), log="y", notch=TRUE,
        main = "relative errors for 'random' eigen-ok 20 x 20 matrix")

## Now  try an example with badly conditioned "random" M matrix...
## ...
## ... (not yet)


### m2ex3() --- The 2x2 example with bad condition , see A3 in ./ex2.R

RE <- re.m2ex3(1e-8, expmList)
sort(RE)# Ward + both sps.H08 are best; s.P.s fair, Eigen (and hybrid): ~1e-9

eps <- 10^-(1:18)
t.m2 <- t(sapply(eps, re.m2ex3, EXPMlist = expmList))
## --> 3 error messages from solve(V), 5 error messages from try(. "R_Eigen" ...)

cbind(sort(apply(log(t.m2),2, median, na.rm=TRUE)))
## 'na.rm=TRUE' needed for Eigen which blows up for the last 3 eps
t.m2.ranks <- sort(rowMeans(apply(t.m2, 1, rank)))
cbind(signif(t.m2.ranks, 3))
## sPs.H08. 2.08
## sPs.H08b 2.08
## Ward     2.25
## s.T.s    5.44
## s.T.sO   5.44
## s.P.s    6.06
## s.P.sO   6.06
## hybrid   7.25
## Eigen    8.33

print(t.m2[, names(t.m2.ranks)[1:8]], digits = 3)
## ==> 1st class: H08 (both) and (but slightly better than)  Ward
##     2nd class  s.T.s and s.P.s
##    "bad" : hybrid and Eigen

if(require(RColorBrewer)) {
    ## Bcol <- brewer.pal(ncol(t.m2),"Dark2")
    Bcol <- brewer.pal(ncol(t.m2),"Set1")
} else {
    ## 7 from Dark2
    ## Bcol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    ##		 "#66A61E", "#E6AB02", "#A6761D")
    ## Rather: those from "Set1"
    Bcol <- c("#E41A1C", "#377EB8", "#4DAF4A",
	      "#984EA3", "#FF7F00", "#FFFF33",
	      "#A65628", "#F781BF", "#999999")
}

matplot(eps, t.m2, type = "b", log = "xy", col=Bcol, lty = 1:9, pch=1:9,
        axes=FALSE, frame = TRUE,
        xlab = expression(epsilon), ylab = "relative error",
        main = expression(expm(A, method == "*") *"  relative errors for  " *
            A == bgroup("[", atop({-1} *"  "* 1, {epsilon^2} *"  "*{-1}), "]")))
legend("bottomright",colnames(t.m2),       col=Bcol, lty = 1:7, pch=1:7,
       inset = 0.02)
if(require("sfsmisc")) {
    sfsmisc::eaxis(1, labels=FALSE)
    sfsmisc::eaxis(1, at = eps[c(TRUE,FALSE)])
    sfsmisc::eaxis(2, labels=FALSE)
    op <- par(las=2)
    sfsmisc::eaxis(2, at = axTicks(2,log=TRUE)[c(TRUE,FALSE,FALSE)])
    par(op)
} else {
    axis(1)
    axis(2)
}

## typical case:
ep <- 1e-10
(me <- m2ex3(ep))
me$expA * exp(1) ## the correct value ; numerically identical to simple matrix:
stopifnot(identical(me$expA * exp(1),
                    rbind(c(  1,  1),
                          c(ep^2, 1))))
## The relative error (matrices):
lapply(expmList, function(EXPM) 1 - EXPM(me$A)/me$expA)

## Average number of correct digits [less "extreme" than plot above]
nDig <- sapply(expmList, function(EXPM) -log10(mean(abs(1 - EXPM(me$A)/me$expA))))
round(nDig, 2)
##   Ward  s.P.s s.P.sO  s.T.s s.T.sO  Eigen hybrid
##  16.26  14.65  14.65  14.65  14.65   6.20   6.39  [AMD Opteron 64-bit]
##    Inf  14.65  14.65  14.65  14.65   6.74   6.33  [Pentium-M (32-bit)]

###--- rnilMat() : random upper triangular (zero-diagonal) nilpotent  n x n matrix

set.seed(17)
m <- rnilMat(10)
as(m, "sparseMatrix")# for nicer printing
E.m <- expm::expm(m, method="Pade")
as(E.m, "sparseMatrix")

(dN <- 9*7*320) # 20160
stopifnot(abs(round(E.m * dN)  -  (E.m * dN)) < 9e-6)
EmN <- matrix(c(dN, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                3*dN, dN, 0, 0, 0, 0, 0, 0, 0, 0,
                352800, 5*dN, dN, 0, 0, 0, 0, 0, 0, 0,
                1018080, 332640, 5*dN, dN, 0, 0, 0, 0, 0, 0,
                2235240, 786240, 292320, 3*dN, dN, 0, 0, 0, 0, 0,
                9368520, 3483480, 1582560, 413280, 181440, dN, 0, 0, 0, 0,
                24676176, 9598680, 5073600, 1562400, 826560, 161280, dN, 0,0,0,
                43730160, 17451000, 10051440, 3430560, 1955520, 504000,
                5*dN, dN, 0, 0,
                68438436, 27747480, 16853760, 6036240, 3638880, 1038240,
                252000, 3*dN, dN, 0,
                119725855, 49165892, 31046760, 11652480, 7198800, 2264640,
                614880, 191520, 3*dN, dN),
              10, 10)

Em.xct <- EmN / dN

stopifnot(all.equal(E.m, Em.xct,
                    check.attributes = FALSE, tol= 1e-13))
re.x <- sapply(expmL.wo.E, function(EXPM) relErr(Em.xct, EXPM(m)))
## with error message from "safe.Eigen"  -->  Eigen is NA here

## result depends quite a bit on platform:
options(digits = 4, width=90)

## Pentium-M 32-bit ubuntu gave
##      Ward     s.P.s    s.P.sO  sPs.H08.  sPs.H08b     s.T.s    s.T.sO    hybrid
## 1.079e-16 4.505e-14 4.503e-14 9.379e-17 9.379e-17 3.716e-17 7.079e-18 1.079e-16
## 32-bit Quad-Core AMD Opteron 2380 (Linux 2.6.30.10-105.2.23.fc11.i686.PAE):
##      Ward     s.P.s    s.P.sO  sPs.H08.  sPs.H08b     s.T.s    s.T.sO    hybrid
## 1.079e-16 4.505e-14 4.503e-14 9.379e-17 9.379e-17 3.716e-17 7.079e-18 1.079e-16

## "Ward77": again more accurate than s+Pade+s, but s+Taylor+s is even more accurate

## but on 64-bit AMD Opterons
##     Ward    s.P.s   s.P.sO sPs.H08. sPs.H08b    s.T.s   s.T.sO   hybrid
## 4.42e-17 3.99e-17 3.99e-17 1.10e-16 1.10e-16 8.44e-17 8.44e-17 4.42e-17
##
## even more astonishing the result on Mac OSX (x86_32_mac; R-forge, R 2.9.0 patch.)
##     Ward    s.P.s   s.P.sO sPs.H08. sPs.H08b    s.T.s   s.T.sO hybrid
## 5.13e-17 3.99e-17 3.99e-17 1.84e-15 1.84e-15 8.44e-17 8.44e-17 5.13e-17

which(is.na(re.x))
(re.x <- re.x[!is.na(re.x)])

stopifnot(re.x[c("Ward", "s.T.s", "s.T.sO")] < 3e-16,
          re.x < 1e-13)# <- 32-bit needed 0.451e-14

cat('Time elapsed: ',(p3 <- proc.time())-p2,'\n') # for ``statistical reasons''
