library(expm)

source(system.file("test-tools.R", package= "expm"), keep.source=FALSE)
## -> assertError()... doExtras
doExtras

## Missing REPROTECT(), till 2014-09-03 [because 'A' is *integer*]:
set.seed(17)
n <- 300
A <- matrix(rbinom(n^2, size=1, prob=0.1), n,n)
A2 <- A %^% 2
for(i in 1:100) {
    A. <- A %^% 2
    if(!isTRUE(all.equal(A2, A.)))
        cat("not equal; i=",i,"\n")
}
## MM: On nb-mm3, I get a different error which shows memory corruption:
##     REAL() can only be applied to a 'numeric', not a 'character'
## or  REAL() can only be applied to a 'numeric', not a 'NULL'

## Check that *large* matrices now work
if(FALSE) ## << even m %^% 2 takes > 20 hours (!!!) [but no longer stops early!]
if(doExtras && require("sfsmisc") &&
   exists("Sys.memGB", "package:sfsmisc", mode="function") &&
   sfsmisc::Sys.memGB() > 50) { ## seems to need   3 x size(m)
    ##
    n <- 46341
    print(as.integer(n^2))# integer overflow
    cat("Creating large matrix 'm' (more than max_int entries):\n ")
    print(system.time(m <- diag(x = (1:n)^3, nrow = n))) # 9.1 sec
    i <- 1:(n-1)
    print(system.time( m[cbind(i,i+1)] <- i ))           # 11.3 sec
    cat("object.size(m): "); print(object.size(m), units="Gb")
    ## 16 Gb (= 17.78 e9 bytes)
    ## This __STILL__ takes hours
    cat("m %^% 2: "); print(system.time(m2 <- m %^% 2))
    ##       user     system    elapsed 
    ## 127199.580   9608.373 137236.405  ==>

    cat("m %^% 4: "); print(system.time(m4 <- m %^% 4)) #
}

