library(expm)

source(system.file("test-tools.R", package= "expm"))## -> assertError()...

## Note that these results are achieved with the default
## settings order=8, method="Pade" -- accuracy could
## presumably be improved still further by some tuning
## of these settings.


## ----------------------------
## Test case 1 from Ward (1977)
## ----------------------------
test1 <- rbind(c(4, 2, 0),
               c(1, 4, 1),
               c(1, 1, 4))
(m1 <- expm(test1, method="Pade"))
(m1O <- expm(test1, method="PadeO"))# very slightly different
(m1T  <- expm(test1,method="Taylor"))
(m1TO <- expm(test1,method="TaylorO"))

## "true" result
m1.t <- matrix(c(147.866622446369, 127.781085523181, 127.781085523182,
                 183.765138646367, 183.765138646366, 163.679601723179,
                 71.797032399996,  91.8825693231832, 111.968106246371), 3,3)
stopifnot(all.equal(m1.t, m1,  check.attrib=FALSE, tol = 1e-13),
          all.equal(m1.t, m1O, check.attrib=FALSE, tol = 1e-13),
          all.equal(m1.t,m1T,  check.attrib=FALSE, tol = 1e-13),
          all.equal(m1.t,m1TO, check.attrib=FALSE, tol = 1e-13),
          all.equal(m1.t, expm(test1,"Ward77"),    tol = 1e-13),
          all.equal(m1.t, expm(test1,"R_Pade"),    tol = 1e-13),
          all.equal(m1.t, expm(test1,"R_Ward77"),  tol = 1e-13))
## -- these agree with ward (1977, p608)
##
m1.2 <- try( expm(test1, "R_Eigen") ) ## 32-bit: gives an error from solve; 64-bit "ok"
if(!inherits(m1.2, "try-error")) {
    stopifnot(all.equal(m1.t, m1.2, check.attrib=FALSE))
    ## but it's less accurate:
    print(all.equal(m1.t, m1.2, check.attrib=FALSE, tol= 1e-12))
    ##-> rel.diff = 6.44e-10
}

##
## ----------------------------
## Test case 2 from Ward (1977)
## ----------------------------
test2 <- t(matrix(c(
                    29.87942128909879, .7815750847907159, -2.289519314033932,
                    .7815750847907159, 25.72656945571064,  8.680737820540137,
                    -2.289519314033932, 8.680737820540137,  34.39400925519054),
                  3, 3))
(m2 <- expm(test2, method="Pade"))
##                   [,1]               [,2]               [,3]
##[1,]   5496313853692357 -18231880972009844 -30475770808580828
##[2,] -18231880972009852  60605228702227024 101291842930256144
##[3,] -30475770808580840 101291842930256144 169294411240859072
## -- which agrees with Ward (1977) to 13 significant figures
(m2O  <- expm(test2, method="PadeO"))
(m2T  <- expm(test2,method="Taylor"))
(m2TO <- expm(test2,method="TaylorO"))

m2.t <- matrix(c(5496313853692216, -18231880972008932, -30475770808579672,
                 -18231880972008928, 60605228702222480, 101291842930249776,
                 -30475770808579672, 101291842930249808, 169294411240850528),
               3, 3)

## -- in this case a very similar degree of accuracy -- even Taylor is ok
stopifnot(all.equal(m2.t, m2, check.attrib=FALSE, tol = 1e-12),
          all.equal(m2.t, m2O,check.attrib=FALSE, tol = 1e-12),
          all.equal(m2.t,m2T, check.attrib=FALSE, tol = 1e-12),
          all.equal(m2.t,m2TO,check.attrib=FALSE, tol = 1e-12),
          all.equal(m2.t, expm(test2,"Ward77"),   tol = 1e-12),
          all.equal(m2.t, expm(test2,"R_Ward77"), tol = 1e-12),
          all.equal(m2.t, expm(test2,"R_Pade"),   tol = 1e-12),
          TRUE)

## ----------------------------
## Test case 3 from Ward (1977)
## ----------------------------
test3 <- t(matrix(c(
    -131, 19, 18,
    -390, 56, 54,
    -387, 57, 52), 3, 3))
(m3 <- expm(test3, method="Pade"))
##                    [,1]                [,2]                [,3]
##[1,] -1.5096441587713636 0.36787943910439874 0.13533528117301735
##[2,] -5.6325707997970271 1.47151775847745725 0.40600584351567010
##[3,] -4.9349383260294299 1.10363831731417195 0.54134112675653534
## -- agrees to 10dp with Ward (1977), p608.
(m3O  <- expm(test3, method="PadeO"))
(m3T  <- expm(test3,method="Taylor"))
(m3TO <- expm(test3,method="TaylorO"))

m3.t <- matrix(c(-1.50964415879218, -5.6325707998812, -4.934938326092,
                 0.367879439109187, 1.47151775849686, 1.10363831732856,
                 0.135335281175235, 0.406005843524598, 0.541341126763207),
               3,3)

stopifnot(all.equal(m3.t, m3,           check.attrib=FALSE, tol = 1e-11),
          all.equal(m3.t, m3T,          check.attrib=FALSE, tol = 1e-11),
          all.equal(m3.t, m3O,          check.attrib=FALSE, tol = 1e-11),
          all.equal(m3.t, m3TO,         check.attrib=FALSE, tol = 1e-11),
          all.equal(m3.t, expm(test3,"R_Eigen"), tol = 1e-11),
          all.equal(m3.t, expm(test3,"Ward77"), tol = 1e-11),
          all.equal(m3.t, expm(test3,"R_Ward"), tol = 1e-11),
          all.equal(m3.t, expm(test3,"R_Pade"), tol = 1e-11),
          TRUE)
## -- in this case, a similar level of agreement with Ward (1977).

## ----------------------------
## Test case 4 from Ward (1977)
## ----------------------------
test4 <-
    array(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-10,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
          dim = c(10, 10))
(m4   <- expm(test4, method="Pade"))
(m4O  <- expm(test4, method="PadeO"))
(m4T  <- expm(test4,method="Taylor"))
(m4TO <- expm(test4,method="TaylorO"))

stopifnot(all.equal(m4  [,10], 1/gamma(10:1), tol=1e-14),
          all.equal(m4O [,10], 1/gamma(10:1), tol=1e-14),
          all.equal(m4T [,10], 1/gamma(10:1), tol=1e-14),
          all.equal(m4TO[,10], 1/gamma(10:1), tol=1e-14),
	  all.equal(m4, m4O, check.attrib=FALSE, tol=5e-15),
	  all.equal(m4, m4T, check.attrib=FALSE, tol=5e-15),
	  all.equal(m4, m4TO,check.attrib=FALSE, tol=5e-15),
          all.equal(m4, expm(test4,"Ward77"), check.attrib=FALSE, tol = 1e-14),
          all.equal(m4, expm(test4,"R_Ward"), check.attrib=FALSE, tol = 1e-14),
          all.equal(m4, expm(test4,"R_Pade"), check.attrib=FALSE, tol = 1e-14),
          max(abs(m4 - expm(test4,"R_Eigen"))) < 1e-7)
## here expm(., EV ) is accurate only to 7 d.p., whereas
##      expm(.,Pade) is correct to at least 14 d.p.
