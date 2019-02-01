## these tests are long, only run if I am me
f <- "~/.R/runalltests"
if(!file.exists(f))
    q(save="no")

## work around for bug in testthat/data.table collaboration
## https://stackoverflow.com/questions/13106018/data-table-error-when-used-through-knitr-gwidgetswww/13131555#13131555
assignInNamespace("cedta.override",
                  c(data.table:::cedta.override,"MFM"),
                  "data.table")

library(testthat)
library(MFM)
library(data.table)
## sharing parameter, largish, so we expect to see a big enough difference to check if working
S <- 100
prior.bin.fn <- function(nT1,s=100,mT1=3) {
    dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
}

makeM <- function(str,usnps) {
    M <- matrix(0,length(str),length(usnps),dimnames=list(NULL,usnps))
    ss <- strsplit(str,"%")
    for(i in seq_along(ss))
        M[i,setdiff(ss[[i]],"0")] <- 1
    return(M)
}
    
if(FALSE) {
    devtools::install_local("~/RP/MTFM")
}

intoverlap <- function(s1,s2) {
    if(length(intersect(s1,s2))>0) {
        1
    } else {
        0
    }
}

threepp <- function(data) {
    ss1 <- strsplit(data$t1.str,"%")
    ss2 <- strsplit(data$t2.str,"%")
    ss3 <- strsplit(data$t3.str,"%")
    data$overlap <- sapply(1:nrow(data), function(i)
        intoverlap(ss1[[i]],ss2[[i]]) +
        intoverlap(ss1[[i]],ss3[[i]]) +
        intoverlap(ss2[[i]],ss3[[i]])) / 3
    data$apr <- prior.bin.fn(data$t1.size) * prior.bin.fn(data$t2.size) *
    prior.bin.fn(data$t3.size) * S^(data$overlap)
    data$apr <- data$apr/sum(data$apr)
    data$pp <- MTFM:::calcpp(data$apr,data$t1t2.logbf)
    data$rpr <- prior.bin.fn(data$t1.size) * prior.bin.fn(data$t2.size) *
    prior.bin.fn(data$t3.size) 
    data$rpr <- data$rpr/sum(data$rpr)
    data$rpp <- MTFM:::calcpp(data$rpr,data$t1t2.logbf)
    return(data)
}
}


################################################################################

context("testing marginalpp")
## load data
testthat::skip_on_travis()
data <- read.table("/scratch/wallace/IL2RA/AD-AC-rep_100-bf-final.txt",header=TRUE,as.is=TRUE)
data <- as.data.table(data)
d1 <- unique(data[,.(t1.logbf,t1.str,t1.size)],by="t1.str")
d2 <- unique(data[,.(t2.logbf,t2.str,t2.size)],by="t2.str")
setnames(d1,c("logbf","str","size"))
setnames(d2,c("logbf","str","size"))
f <- function(d) {
    d$pr <- prior.bin.fn(d$size)
    d$pp <- MTFM:::calcpp(d$pr,d$logbf)
    d
}
d1 <- f(d1)
d2 <- f(d2)

## time 3 way combinations
nrow(d1)
nrow(d2)
library(microbenchmark)

load_all()
microbenchmark( marginalpp(STR=list(d1$str,d2$str,d2$str),
                     ABF=list(d1$logbf,d2$logbf,d2$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr),
                     PP=list(d1$pp,d2$pp,d2$pp),
                     kappa=S,
                     p0=prior.bin.fn(0)))
##      min       lq     mean   median       uq      max neval
## 91.38261 93.72761 98.81083 94.58032 97.14498 394.4529   100
 ##     min       lq     mean   median       uq      max neval // store pp1 * pp2
## 84.45908 87.26599 90.25199 88.13486 91.38183 121.6457   100
 ##      min       lq     mean   median       uq    max neval // do idx1 -1 at beginning
 ## 82.46985 85.51192 88.31313 86.53107 89.59482 121.42   100

microbenchmark( marginalpp(STR=list(d1$str,d2$str,d2$str),
                     ABF=list(d1$logbf,d2$logbf,d2$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr),
                     PP=list(d1$pp,d2$pp,d2$pp),
                     kappa=S,
                     p0=prior.bin.fn(0)),
               marginallogpp(STR=list(d1$str,d2$str,d2$str),
                     ABF=list(d1$logbf,d2$logbf,d2$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr),
                     PP=list(log(d1$pp),log(d2$pp),log(d2$pp)),
                     kappa=S,
                     p0=prior.bin.fn(0)))

## time 4 way combinationsmarginalpp(STR=list(d1$str,d2$str,d2$str),
                     r1 <- marginalpp(STR=list(d1$str,d2$str,d2$str),
 ABF=list(d1$logbf,d2$logbf,d2$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr),
                     PP=list(d1$pp,d2$pp,d2$pp),
                     kappa=S,
                     p0=prior.bin.fn(0))
               r2 <- marginallogpp(STR=list(d1$str,d2$str,d2$str),
                     ABF=list(d1$logbf,d2$logbf,d2$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr),
                     PP=list(log(d1$pp),log(d2$pp),log(d2$pp)),
                     kappa=S,
                     p0=prior.bin.fn(0))
lapply(r1$shared.pp,summary)
lapply(r2$shared.pp,summary)

microbenchmark( marginalpp(STR=list(d1$str,d2$str,d2$str,d1$str),
                     ABF=list(d1$logbf,d2$logbf,d2$logbf,d1$logbf),
                     pr=list(d1$pr,d2$pr,d2$pr,d1$pr),
                     PP=list(d1$pp,d2$pp,d2$pp,d1$pp),
                     kappa=S,
                     p0=prior.bin.fn(0)),
               times=10)
 ##      min       lq     mean   median       uq      max neval
 ## 9.289625 9.358305 9.473313 9.456111 9.589571 9.705593    10
 ##      min      lq     mean   median       uq      max neval
## 8.367299 8.59717 8.664065 8.644912 8.741663 8.916749    10
 ##   min       lq     mean   median       uq      max neval
 ## 7.681974 7.709554 7.880879 7.925056 7.961007 8.137236    10
################################################################################


## time stroverlap
load_all()
microbenchmark(MTFM:::stroverlap(c(1,3,8,10),c(2,6,10)),
               MTFM:::strint(c(1,3,8,10),c(2,6,10)),
               MTFM:::stroverlap(c(1,3,8,10),c(1,6,10)),
               MTFM:::strint(c(1,3,8,10),c(1,6,10)) ) 

## make 4 way combinations
    data <- expand.grid(1:nrow(marg$d1),1:nrow(marg$d1),1:nrow(marg$d2),1:nrow(marg$d2))
    data <- cbind(data,marg$d1[adata$Var1,.(str,logbf,size)])
    data <- cbind(data,marg$d1[adata$Var2,.(str,logbf,size)])
    data <- cbind(data,marg$d2[adata$Var3,.(str,logbf,size)])
    data <- cbind(data,marg$d2[adata$Var4,.(str,logbf,size)])
dim(data)
colnames(data)[-c(1:4)] <- c("t1.str","t1.logbf","t1.size",
                                  "t2.str","t2.logbf","t2.size",
                                  "t3.str","t3.logbf","t3.size","t4.str","t4.logbf","t4.size" )
    data$t1t2.logbf <- data$t1.logbf + data$t2.logbf + data$t3.logbf- 46.6
head(data)
    data <- threepp(data)


dt <- as.data.table(data)
d1 <- (dt[,.(t1.logbf=unique(t1.logbf),t1.size=unique(t1.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
          by=c("t1.str")])
d2 <- (dt[,.(t2.logbf=unique(t2.logbf),t2.size=unique(t2.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
          by=c("t2.str")])
d3 <- (dt[,.(t3.logbf=unique(t3.logbf),t3.size=unique(t3.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
          by=c("t3.str")])
setnames(d1,c("str","logbf","size","spp.orig","mpp.orig"))
    setnames(d2,c("str","logbf","size","spp.orig","mpp.orig"))
    setnames(d3,c("str","logbf","size","spp.orig","mpp.orig"))
    d1 <- unique(d1,by="str")
    d2 <- unique(d2,by="str")
    d3 <- unique(d3,by="str")
    d1 <- d1[str!="0",]
    d2 <- d2[str!="0",]
    d3 <- d3[str!="0",]
    d1$pr <- prior.bin.fn(d1$size)
    d2$pr <- prior.bin.fn(d2$size)
    d3$pr <- prior.bin.fn(d3$size)
    d1$mpp.direct <- MTFM:::calcpp(d1$pr,d1$logbf)
    d2$mpp.direct <- MTFM:::calcpp(d2$pr,d2$logbf)
    d3$mpp.direct <- MTFM:::calcpp(d3$pr,d3$logbf)
    ss1 <- strsplit(data$t1.str,"%")
    ss2 <- strsplit(data$t2.str,"%")
    ss3 <- strsplit(data$t3.str,"%")
    usnps <- setdiff(unique(unlist(c(ss1,ss2,ss3))),"0")
    length(usnps)
    M1 <- makeM(d1$str,usnps)
    M2 <- makeM(d2$str,usnps)
    M3 <- makeM(d3$str,usnps)
    PPS <- marginalpp(STR=list(d1$str,d2$str,d3$str),
                     ABF=list(d1$logbf,d2$logbf,d3$logbf),
                     pr=list(d1$pr,d2$pr,d3$pr),
                     PP=list(d1$mpp.direct,d2$mpp.direct,d3$mpp.direct),
                     kappa=S,
                     p0=prior.bin.fn(0))
microbenchmark::microbenchmark( marginalpp(STR=list(d1$str,d2$str,d3$str),
                     ABF=list(d1$logbf,d2$logbf,d3$logbf),
                     pr=list(d1$pr,d2$pr,d3$pr),
                     PP=list(d1$mpp.direct,d2$mpp.direct,d3$mpp.direct),
                     kappa=S,
                     p0=prior.bin.fn(0)))




## now perform same analysis as above
    amarg <- margthree(adata)

test_that("marginalpp works for 3 diseases", {
    testthat::skip_on_cran()
    expect_equal(amarg$d1$spp.orig,amarg$d1$spp.models)
    expect_equal(amarg$d1$spp.orig,amarg$d1$spp.str)
})


