#### Examples where we know the result "exactly"

library(expm)

options(digits = 4, width = 99, keep.source = FALSE)

mSource <- function(file, ...)
    source(system.file(file, ..., package = "expm", mustWork=TRUE))
mSource("test-tools.R")## -> assertError(), rMat(), .. doExtras
mSource("demo", "exact-fn.R")# -> nilA3(), facMat(), m2ex3(), ....
doExtras

re.nilA3 <- function(xyz, EXPMlist)
{
    stopifnot(is.list(EXPMlist))
    r <- do.call(nilA3, as.list(xyz))
    sapply(EXPMlist, function(Efn) relErr(r$expA, Efn(r$A)))
}

re.facMat <- function(n, EXPMlist, rFUN = rnorm, ...)
{
    stopifnot(is.list(EXPMlist))
    r <- facMat(n, rFUN, ...)
    vapply(EXPMlist, function(EXPM) {
	ct <- system.time(E <- EXPM(r$A), gcFirst = FALSE)[[1]]
	c(relErr = relErr(r$expA, E), c.time = ct)
    }, double(2))
}

re.m2ex3 <- function(eps, EXPMlist)
{
    stopifnot(is.list(EXPMlist))
    r <- m2ex3(eps)
    sapply(EXPMlist, function(EXPM) relErr(r$expA, EXPM(r$A)))
}

## check 1x1 matrices
stopifnot(
    ## these had failed before 2017-03-28 (= Liselotte's 58-th birthday):
    all.equal(as.matrix(sqrtm(matrix(4))),  matrix(2)) ,
    all.equal(as.matrix(logm (matrix(pi))), matrix(log(pi))) ,
  ## these had "always" worked :
    all.equal(as.matrix(expm (matrix(0))), matrix(1)) ,
    all.equal(as.matrix(expm (matrix(1))), matrix(exp(1)))
)


set.seed(321)
re <- replicate(1000,
                c(re.nilA3(rlnorm(3),list(function(x)expm(x,"Pade"))),
                  re.nilA3(rnorm(3), list(function(x)expm(x,"Pade")))))

summary(t(re))
stopifnot(rowMeans(re) < 1e-15,
          apply(re, 1, quantile, 0.80) < 1e-16,
          apply(re, 1, quantile, 0.90) < 2e-15,
          apply(re, 1, max) < c(4e-14, 6e-15))

showProc.time()

## Check *many* random nilpotent  3 x 3  matrices:
set.seed(321)
RE <- replicate(1000,
                c(re.nilA3(rlnorm(3), list(function(x) expm(x, "Ward77"))),
                  re.nilA3(rnorm(3),  list(function(x) expm(x, "Ward77")))))
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
## nb-mm3; ditto lynne (both x64), 2014-09-11:
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## -0.8422  0.0000  0.0000  0.0725  0.1067  1.2205
## 32-bit [i686, florence, Linx 3.14.8-100.fc19..]:
##  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
## -0.62    0.00    0.00    1.36    0.14   95.93


showProc.time()

###--- A second group --- where we know the diagonalization of A ---

if(!require("Matrix"))
    q('no')
##  ------  the rest really uses 'Matrix'
##---> now use  expm::expm()  since Matrix has its own may mask the expm one
##              ^^^^^^^^^^^^

(osV <- abbreviate(gsub("[^[:alnum:]]", '', sub("\\(.*", '', osVersion)), 12))
if(!dev.interactive(TRUE)) pdf(paste0("expm_exact-ex_", osV, ".pdf"), width = 9, height=5)
##
myRversion <- paste(R.version.string, "--", osVersion)
if((mach <- Sys.info()[["machine"]]) != "x86_64")
    myRversion <- paste0(myRversion, "_", mach)
if(!capabilities("long.double"))
    myRversion <- paste0(myRversion, "_no_LDbl")
myRversion
## in plots use   mtext(myRversion, adj=1, cex=3/4)


## rMat() relies on Matrix::rcond():
## Now with the change default rcondMin, this "works"
R40 <- rMat(40)
R80 <- rMat(80)
showProc.time()

expm.safe.Eigen <- function(x, silent = FALSE) {
    r <- try(expm::expm(x, "R_Eigen"), silent = silent)
    if(inherits(r, "try-error")) NA else r
}

## the S4 generic
Matrix::expm
## the dgeMatrix method - adapted to Matrix changes, had *more versatile*  ..2dge() :
expm.Matr.dge <- function(x) getDataPart(getMethod("expm", "dgeMatrix"))(Matrix::.m2dense(x))

expmList <-
    list(Matr = Matrix::expm,
         Matr.d   = expm.Matr.dge,
         Ward     = function(x) expm::expm(x, "Ward77"),
         Ward1b   = function(x) expm::expm(x, "Ward77", preconditioning = "1bal"),
         RWard6   = function(x) expm::expm(x, "R_Ward77", order = 6),
         RWard7   = function(x) expm::expm(x, "R_Ward77", order = 7),
         RWard8   = function(x) expm::expm(x, "R_Ward77", order = 8), # default
         RWard9   = function(x) expm::expm(x, "R_Ward77", order = 9),
	 s.P.s7   = function(x) expm::expm(x, "Pade", order = 7),
	 s.POs7   = function(x) expm::expm(x, "PadeO",order = 7),
	 s.P.s8   = function(x) expm::expm(x, "Pade" ,order = 8), # default
	 s.POs8   = function(x) expm::expm(x, "PadeO",order = 8), # default
	 s.P.s9   = function(x) expm::expm(x, "Pade", order = 9),
	 s.POs9   = function(x) expm::expm(x, "PadeO",order = 9),
	 s.P.sRBS = function(x) expm::expm(x, "PadeRBS"),
	 Rs.P.s7  = function(x) expm::expm(x, "R_Pade", order = 7),
	 Rs.P.s8  = function(x) expm::expm(x, "R_Pade", order = 8), # default
	 Rs.P.s9  = function(x) expm::expm(x, "R_Pade", order = 9),
	 sPs.H08. = function(x) expm:: expm.Higham08(x, balancing=FALSE),
	 sPs.H08b = function(x) expm:: expm.Higham08(x, balancing= TRUE),
	 ## AmHi09.06= function(x) expm:::expm.AlMoHi09(x, p =  6),
	 AmHi09.07= function(x) expm:::expm.AlMoHi09(x, p =  7),
	 AmHi09.08= function(x) expm:::expm.AlMoHi09(x, p =  8), # default -- really sub optimal
	 AmHi09.09= function(x) expm:::expm.AlMoHi09(x, p =  9),
	 AmHi09.10= function(x) expm:::expm.AlMoHi09(x, p = 10),
	 AmHi09.11= function(x) expm:::expm.AlMoHi09(x, p = 11),
	 AmHi09.12= function(x) expm:::expm.AlMoHi09(x, p = 12),
	 AmHi09.13= function(x) expm:::expm.AlMoHi09(x, p = 13),
	 s.T.s7   = function(x) expm::expm(x, "Taylor", order = 7),
	 s.TOs7   = function(x) expm::expm(x, "TaylorO",order = 7),
	 s.T.s8   = function(x) expm::expm(x, "Taylor", order = 8), # default
	 s.TOs8   = function(x) expm::expm(x, "TaylorO",order = 8), # default
	 s.T.s9   = function(x) expm::expm(x, "Taylor", order = 9),
	 s.TOs9   = function(x) expm::expm(x, "TaylorO",order = 9),
	 Eigen    = expm.safe.Eigen, # "R_Eigen"
	 hybrid   = function(x) expm::expm(x, "hybrid")
	 )


## set.seed(12)
## facMchk <- replicate(if(doExtras) 100 else 20, facMat(20, rnorm))
set.seed(12)
fRE <- replicate(if(doExtras) 100 else 20,
                 re.facMat(20, expmList)) # if(doExtras) gives one "No Matrix found ..." warning
nDig <- -log10(t(fRE["relErr",,]))
cat("Number of correct decimal digits for facMat(20, rnorm):\n")
t(apply(nDig, 2, summary))

## Now look at that:
eaxis <- if(requireNamespace("sfsmisc")) sfsmisc::eaxis else axis
op <- par(mar=.1 + c(5,4 + 1.5, 4,2))
boxplot(t(fRE["relErr",,]), log="x", xaxt="n",
        notch=TRUE, ylim = c(8e-16, 4e-9),
        horizontal=TRUE, las = 1,
        main = "relative errors for 'random' eigen-ok 20 x 20 matrix")
eaxis(1); grid(lty = 3);  mtext(myRversion, adj=1, cex=3/4)
par(op)

showProc.time()

if(doExtras) withAutoprint({ # also "large" n = 100 ------------------------------------------
str(rf100 <- replicate(20, re.facMat(100, expmList)))
1000*t(apply(rf100["c.time",,], 1, summary))
## lynne {Linux 2.6.34.7-56.fc13.x86_64 --- AMD Phenom II X4 925}:
##          Min. 1st Qu. Median  Mean 3rd Qu. Max.
## Ward       23      24   24.5  24.4    25.0   25
## s.P.s     107     109  109.0 109.0   109.0  112
## s.POs    188     190  191.0 192.0   193.0  198
## s.P.sRBS   17      18   19.0  18.9    19.2   21
## sPs.H08.   15      17   18.0  17.6    18.0   19
## sPs.H08b   18      18   19.0  23.4    20.0  107
## s.T.s      44      45   45.0  45.6    46.0   48
## s.TOs     96      98   99.0 100.0   100.0  116
## Eigen      18      19   20.0  24.4    21.0  109
## hybrid     40      42   42.0  47.1    44.0  133

nDig <- -log10(t(rf100["relErr",,]))
cat("Number of correct decimal digits for facMat(100, rnorm):\n")
t(apply(nDig, 2, summary))

##--> take out the real slow ones for the subsequent tests:
(not.slow <- grep("^s\\.[PT]", names(expmList), invert = TRUE, value = TRUE))

set.seed(18) ## 12 replicates is too small .. but then it's too slow otherwise:
rf400 <- replicate(12, re.facMat(400, expmList[not.slow]))

showProc.time()
1000*t(apply(rf400["c.time",,], 1, summary))
## lynne:
##          Min. 1st Qu. Median Mean 3rd Qu. Max.
## Ward     1740    1790   1830 1820    1860 1900
## s.P.sRBS 1350    1420   1440 1430    1450 1460
## sPs.H08. 1020    1030   1130 1140    1210 1290
## sPs.H08b 1120    1130   1220 1220    1300 1390
## Eigen     962     977    989  992    1000 1030
## hybrid   2740    2800   2840 2840    2890 2910

nDig <- -log10(t(rf400["relErr",,]))
cat("Number of correct decimal digits for (12 rep. of) facMat(400, rnorm):\n")
t(apply(nDig, 2, summary))

}) else { # *not* (doExtras) -----------------------------------------------------------------
    ## times (in millisec):
    print(1000*t(apply(fRE["c.time",,], 1, summary)))
}

## Now  try an example with badly conditioned "random" M matrix...
## ...
## ...   (not yet -- TODO?)

### m2ex3() --- The 2x2 example with bad condition , see A3 in ./ex2.R

RE <- re.m2ex3(1e-8, expmList)
sort(RE)# Ward + both sps.H08 are best; s.P.s fair, Eigen (and hybrid): ~1e-9

eps <- 10^-(1:18)
t.m2 <- t(sapply(eps, re.m2ex3, EXPMlist = expmList))
## --> 3 error messages from solve(V), 5 error messages from try(. "R_Eigen" ...)

showProc.time()

cbind(sort(apply(log(t.m2),2, median, na.rm=TRUE)))
## 'na.rm=TRUE' needed for Eigen which blows up for the last 3 eps
t.m2.ranks <- sort(rowMeans(apply(t.m2, 1, rank)))
cbind(signif(t.m2.ranks, 3))
## lynne (x86_64, Linux 3.14.4-100; Intel i7-4765T), 2014-09:
## sPs.H08.   2.67
## sPs.H08b   2.67
## s.P.sRBS   3.06
## Ward       4.03
## AmHi09.13  4.33 <<- still not close to H08 !
## AmHi09.12  5.86
## s.T.s      8.33
## s.TOs     8.33
## s.P.s      9.11
## s.POs     9.11
## hybrid    10.80
## AmHi09.10 11.70 << astonishingly bad
## Eigen     12.60
## AmHi09.08 13.10
## AmHi09.06 14.40

print(t.m2[, names(t.m2.ranks)[1:8]], digits = 3)
## ==> 1st class: H08 (both) and (but slightly better than)  Ward
##     2nd class  s.T.s and s.P.s
##    "bad" : hybrid and Eigen

## ??? AmHi09 - methods, up to order = 10 are worse !
if(require(RColorBrewer)) {
    ## Bcol <- brewer.pal(ncol(t.m2),"Dark2")
    Bcol <- brewer.pal(min(9, ncol(t.m2)), "Set1")
    Bcol <- Bcol[sqrt(colSums(col2rgb(Bcol)^2)) < 340]
    ## FIXME: more colors ==> ~/R/MM/GRAPHICS/color-palettes.R
} else {
    ## 7 from Dark2
    ## Bcol <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    ##		 "#66A61E", "#E6AB02", "#A6761D")
    ## Rather: those from "Set1"
    Bcol <- c("#E41A1C", "#377EB8", "#4DAF4A",
	      "#984EA3", "#FF7F00", # too bright: "#FFFF33",
	      "#A65628", "#F781BF", "#999999")
}

matplot(eps, t.m2, type = "b", log = "xy", col=Bcol, lty = 1:9, pch=1:15,
        axes=FALSE, frame = TRUE,
        xlab = expression(epsilon), ylab = "relative error",
        main = expression(expm(A, method == "*") *"  relative errors for  " *
            A == bgroup("[", atop({-1} *"  "* 1, {epsilon^2} *"  "*{-1}), "]")))
legend("bottomright",colnames(t.m2),       col=Bcol, lty = 1:9, pch=1:15,
       inset = 0.02)
if(requireNamespace("sfsmisc")) {
    sfsmisc::eaxis(1, labels=FALSE)
    sfsmisc::eaxis(1, at = eps[c(TRUE,FALSE)])
    sfsmisc::eaxis(2)
    ## sfsmisc::eaxis(2, labels=FALSE)
    ## op <- par(las=2)
    ## sfsmisc::eaxis(2, at = axTicks(2,log=TRUE)[c(TRUE,FALSE,FALSE)])
    ## par(op)
} else {
    axis(1)
    axis(2)
}

## typical case:
ep <- 1e-10
(me <- m2ex3(ep))
me$expA * exp(1) ## the correct value ; numerically identical to simple matrix:
## identical() not fulfilled e.g. on Solaris
stopifnot(all.equal(me$expA * exp(1),
		    rbind(c(  1,  1),
			  c(ep^2, 1)),
		    tolerance = 1e-14))
## The relative error (matrices):
lapply(expmList, function(EXPM) 1 - EXPM(me$A)/me$expA)

## Average number of correct digits [less "extreme" than plot above]
nDig <- sapply(expmList, function(EXPM) -log10(mean(abs(1 - EXPM(me$A)/me$expA))))
round(nDig, 2)
##   Ward  s.P.s s.POs  s.T.s s.TOs  Eigen hybrid
##  16.26  14.65  14.65  14.65  14.65   6.20   6.39  [AMD Opteron 64-bit]
##    Inf  14.65  14.65  14.65  14.65   6.74   6.33  [Pentium-M (32-bit)]

showProc.time()

###--- rnilMat() : random upper triangular (zero-diagonal) nilpotent  n x n matrix

set.seed(17)
m <- rnilMat(10)
(m. <- as(m,"sparseMatrix"))# for nicer printing - and more *below*
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
          all.equal(E.m, Em.xct, check.attributes = FALSE, tolerance = 0)
stopifnot(all.equal(E.m, Em.xct, check.attributes = FALSE, tolerance= 1e-13))

re.x <- sapply(expmList, function(EXPM) relErr(Em.xct, EXPM(m)))
## with error message from "safe.Eigen"  -->  Eigen is NA here
## result depends quite a bit on platform
which(is.na(re.x)) # gives  c(Eigen = 16L)  (but not everywhere ?!?)
(re.x <- re.x[!is.na(re.x)])

## Pentium-M 32-bit ubuntu gave
##      Ward     s.P.s    s.POs  sPs.H08.  sPs.H08b     s.T.s    s.TOs    hybrid
## 1.079e-16 4.505e-14 4.503e-14 9.379e-17 9.379e-17 3.716e-17 7.079e-18 1.079e-16
## 32-bit Quad-Core AMD Opteron 2380 (Linux 2.6.30.10-105.2.23.fc11.i686.PAE):
##      Ward     s.P.s    s.POs  sPs.H08.  sPs.H08b     s.T.s    s.TOs    hybrid
## 1.079e-16 4.505e-14 4.503e-14 9.379e-17 9.379e-17 3.716e-17 7.079e-18 1.079e-16

## "Ward77": again more accurate than s+Pade+s, but s+Taylor+s is even more accurate

## but on 64-bit AMD Opterons
##     Ward    s.P.s   s.POs sPs.H08. sPs.H08b    s.T.s   s.TOs   hybrid
## 4.42e-17 3.99e-17 3.99e-17 1.10e-16 1.10e-16 8.44e-17 8.44e-17 4.42e-17
##
## even more astonishing the result on Mac OSX (x86_32_mac; R-forge, R 2.9.0 patch.)
##     Ward    s.P.s   s.POs sPs.H08. sPs.H08b    s.T.s   s.TOs hybrid
## 5.13e-17 3.99e-17 3.99e-17 1.84e-15 1.84e-15 8.44e-17 8.44e-17 5.13e-17

## 2014-09: AmHi09 are very good (64bit: 8e-17) for p >= 8 (p=6 has 1.5e-11)

not.09.06 <- which(names(re.x) != "AmHi09.06")
stopifnot(re.x[c("Ward", "s.T.s8", "s.TOs8")] < 3e-16,
          ## re.x[["AmHi09.06"]] < 9e-11, # x64 & 686(lnx): = 1.509e-11
          re.x[not.09.06] < 4e-13)# max: 686(32b): 4.52e-14, x64(lnx): 1.103e-16

##-- Looking at *sparse* matrices: [C,Fortran "dense" code based methods will fail]:
(meths <- eval(formals(expm)$method))
ems <- sapply(meths, function(met)
              tryCatch(expm::expm(m., method=met), error=identity))
ok <- !sapply(ems, is, class2="error")
meths[ok] # now most... are

showProc.time()

## Complex Matrices
re.facMat.Z <- function(n, EXPMlist, rFUN = function(n) rnorm(n) + 1i*rnorm(n), ...)
{
    stopifnot(is.list(EXPMlist))
    r <- facMat(n, rFUN, ...)
    vapply(EXPMlist, function(EXPM) {
	ct <- system.time(E <- EXPM(r$A), gcFirst = FALSE)[[1]]
	c(relErr = relErr(r$expA, E), c.time = ct)
    }, double(2))
}

(nmL <- names(expmList))
##  [1] "Matr"      "Matr.d"    "Ward"      "Ward1b"    "s.P.s"     "s.POs"     "s.P.s7"
##  [8] "s.POs7"    "s.P.s9"    "s.POs9"    "s.P.sRBS"  "sPs.H08."  "sPs.H08b"  "AmHi09.06"
## [15] "AmHi09.07" "AmHi09.08" "AmHi09.09" "AmHi09.10" "AmHi09.12" "AmHi09.13" "s.T.s"
## [22] "s.TOs"     "s.T.s7"    "s.TOs7"    "s.T.s9"    "s.TOs9"    "Eigen"     "hybrid"

## dropping "Matr", "Matr.d" (which gives "dgeMatrix" currently --> mean(.) fails ...
##          "Ward" "Ward1b" and "hybrid"  error "not a numeric Matrix"
##          "AmHi09": C code currently only for double precision  ((FIXME))
(cplxN <- grep("^(Matr|Ward|hybr|AmHi09|s\\.[PT])", nmL, invert = TRUE, value = TRUE))

rr <- re.facMat.Z(4, expmList[cplxN])

set.seed(47)
fREZ <- replicate(if(doExtras) 64 else 12,
                 re.facMat.Z(15, expmList[cplxN]))
nDig <- -log10(t(fREZ["relErr",,]))
cat("Number of correct decimal digits for facMat(20, rnorm + i*rnorm):\n")
t(apply(nDig, 2, summary))

## times (in millisec):
print(1000*t(apply(fREZ["c.time",,], 1, summary)))

## Now look at that:
op <- par(mar=.1 + c(5,4 + 1.5, 4,2))
boxplot(t(fREZ["relErr",,]), log="x", xaxt="n",
        notch=TRUE, # ylim = c(8e-16, 4e-9),
        horizontal=TRUE, las = 1,
        main = "relative errors for 'random' <complex> eigen-ok 20 x 20 matrix")
eaxis(1); grid(lty = 3);  mtext(myRversion, adj=1, cex=3/4)
par(op)



showProc.time()
