#### Example matrices from the Matlab demos

library(expm)

source(system.file("test-tools.R", package= "expm"))## -> assertError()...

## --- 1 ---
## Here, all three {eigen; Taylor; Pade(Scaling & Squaring)} should do well

A1 <- rbind(0:2,
            c(0.5, 0, 1),
            2:0)
A1
ml1 <- lapply(c(4:10,20),
              function(order) expm(A1, "Pade", order=order))
for(k in seq_len(length(ml1) - 1))
    stopifnot(all.equal(ml1[[k]], ml1[[k + 1]], tol = 1e-12))

for(k in seq_len(length(ml1) - 1)) {
    print(all.equal(ml1[[k]], ml1[[k + 1]], tol = 0))
}

mA1 <- ml1[[4]]
stopifnot(all.equal(mA1,
                    matrix(c(5.3090812852106, 2.8087900904073, 5.1737460019740,
                             4.0012030182399, 2.8845155413486, 4.0012030182399,
                             5.5778402926177, 3.1930144369526, 5.7131755758543),
                           3, 3),
                           check.attrib = FALSE, tol = 1e-11))


## --- 2 ---
## Here, Taylor "fails":

## A matrix where the terms in the Taylor series become very large
## before they go to zero.
A2 <- rbind(c(-147, 72),
            c(-192, 93))
A2
(mA2 <- expm(A2, method="Pade"))
stopifnot(all.equal(mA2,
                    matrix(c(-0.099574136735459, -0.199148273470915,
                              0.074680602551593 , 0.149361205103183),
                           2, 2), check.attrib = FALSE, tol = 1e-11))
mA2.T <- expm(A2, method = "Taylor")
stopifnot(all.equal(mA2, mA2.T, tol=1e-10))
all.equal(mA2, mA2.T, tol=1e-14)#-> 3.2e-12  {MM: I think that's pretty good}

## --- 3 ---
## Here, Eigenvalues  must fail ("not a full set"):
A3 <- rbind(c(-1, 1),
            c(0, -1))
(mA3 <- expm(A3, method="Pade"))
assertError(expm(mA3, method="R_Eigen"))
em1 <- exp(-1)
stopifnot(all.equal(mA3, ## and the exact solution:
                    matrix(c(em1, 0, em1, em1), 2, 2),
                    check.attrib = FALSE, tol = 1e-14))

## using 'eps' instead of 0 :
## ---> see m2ex3() etc in ./exact-ex.R


## --- 4 ---
## Here, some version of do_expm() failed:
(m <- matrix(c(0,2:0),2))
## Eigenvalue decomposition:
d <- c(sqrt(2), -sqrt(2))
V <- rbind(c(sqrt(1/3), -sqrt(1/3)),
           c(sqrt(2/3),  sqrt(2/3)))
## ==>
IV <- rbind(c( sqrt(3/4), sqrt(3/8)),
            c(-sqrt(3/4), sqrt(3/8)))
stopifnot(all.equal(V %*% IV, diag(2)))
em.true <- V %*% (exp(d) * IV)
stopifnot(all.equal(em.true, expm::expm(m)),
          all.equal(em.true, expm::expm(m,"Pade"), check.attrib=FALSE))
