## these tests are long, only run if I am me
f <- "~/.R/runalltests"
if(!file.exists(f))
    q(save="no")

## work around for bug in testthat/data.table collaboration
## https://stackoverflow.com/questions/13106018/data-table-error-when-used-through-knitr-gwidgetswww/13131555#13131555
assignInNamespace("cedta.override",
                  c(data.table:::cedta.override,"MTFM"),
                  "data.table")

library(testthat)
library(MTFM)
library(data.table)
## sharing parameter, largish, so we expect to see a big enough difference to check if working
S <- 100
prior.bin.fn <- function(nT1,s=100,mT1=3) {
    dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
}


bothpp <- function(data) {
    ss1 <- strsplit(data$t1.str,"%")
    ss2 <- strsplit(data$t2.str,"%")
    data$overlap <- sapply(1:nrow(data), function(i)
        length(intersect(ss1[[i]],ss2[[i]]))>0)
    data$apr <- prior.bin.fn(data$t1.size) * prior.bin.fn(data$t2.size) * ifelse(data$overlap,S,1)
    data$apr <- data$apr/sum(data$apr)
    data$pp <- MTFM:::calcpp(data$apr,data$t1t2.logbf)
    data$rpr <- prior.bin.fn(data$t1.size) * prior.bin.fn(data$t2.size)
    data$rpr <- data$rpr/sum(data$rpr)
    data$rpp <- MTFM:::calcpp(data$rpr,data$t1t2.logbf)
    return(data)
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

margpp <- function(data) {
    dt <- as.data.table(data)
    d1 <- dt[,.(t1.logbf=unique(t1.logbf),t1.size=unique(t1.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
              by=c("t1.str")]
    d2 <- dt[,.(t2.logbf=unique(t2.logbf),t2.size=unique(t2.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
              by=c("t2.str")]
    setnames(d1,c("str","logbf","size","spp.orig","mpp.orig"))
    setnames(d2,c("str","logbf","size","spp.orig","mpp.orig"))
    d1 <- unique(d1,by="str")
    d2 <- unique(d2,by="str")
    d1 <- d1[str!="0",]
    d2 <- d2[str!="0",]
    d1$pr <- prior.bin.fn(d1$size)
    d2$pr <- prior.bin.fn(d2$size)
    d1$mpp.direct <- MTFM:::calcpp(d1$pr,d1$logbf)
    d2$mpp.direct <- MTFM:::calcpp(d2$pr,d2$logbf)
    ss1 <- strsplit(data$t1.str,"%")
    ss2 <- strsplit(data$t2.str,"%")
    usnps <- unique(unlist(c(ss1,ss2)))
    length(usnps)
    M1 <- makeM(d1$str,usnps)
    M2 <- makeM(d2$str,usnps)
    PP <- marginalpp2(M1,M2,d1$logbf,d2$logbf,d1$pr,d2$pr,S,prior.bin.fn(0))
    d1$spp.new <- as.numeric(PP$shared.pp1)[-1]
    d1$mpp.new <- as.numeric(PP$single.pp1)[-1]
    return(list(d1=d1,d2=d2))
}


margpp2 <- function(data) {
    dt <- as.data.table(data)
    d1 <- (dt[,.(t1.logbf=unique(t1.logbf),t1.size=unique(t1.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
              by=c("t1.str")])
    d2 <- (dt[,.(t2.logbf=unique(t2.logbf),t2.size=unique(t2.size),pp.orig=sum(pp),rpp.orig=sum(rpp)),
              by=c("t2.str")])
    setnames(d1,c("str","logbf","size","spp.orig","mpp.orig"))
    setnames(d2,c("str","logbf","size","spp.orig","mpp.orig"))
    d1 <- unique(d1,by="str")
    d2 <- unique(d2,by="str")
    d1$pr <- prior.bin.fn(d1$size)
    d2$pr <- prior.bin.fn(d2$size)
    d1$mpp.direct <- MTFM:::calcpp(d1$pr,d1$logbf)
    d2$mpp.direct <- MTFM:::calcpp(d2$pr,d2$logbf)

    ss1 <- strsplit(data$t1.str,"%")
    ss2 <- strsplit(data$t2.str,"%")
    usnps <- unique(unlist(c(ss1,ss2)))
    length(usnps)
    M1 <- makeM(d1$str,usnps)
    M2 <- makeM(d2$str,usnps)

    PPM <- marginalpp.models(list(M1,M2),list(d1$logbf,d2$logbf),list(d1$pr,d2$pr),S,prior.bin.fn(0))
    PPS <- marginalpp(list(d1$str,d2$str),list(d1$logbf,d2$logbf),list(d1$pr,d2$pr),S,prior.bin.fn(0))
    d1$spp.models <- as.numeric(PPM$shared.pp[[1]])[-1]
    d1$mpp.models <- as.numeric(PPM$single.pp[[1]])[-1]
    d1$spp.str <- as.numeric(PPS$shared.pp[[1]])[-1]
    d1$mpp.str <- as.numeric(PPS$single.pp[[1]])[-1]
    return(list(d1=d1,d2=d2))
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
margthree <- function(data) {
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
    PPM <- marginalpp.models(M=list(M1,M2,M3),
                     ABF=list(d1$logbf,d2$logbf,d3$logbf),
                     pr=list(d1$pr,d2$pr,d3$pr),
                     kappa=S,
                     p0=prior.bin.fn(0))
    PPS <- marginalpp(STR=list(d1$str,d2$str,d3$str),
                     ABF=list(d1$logbf,d2$logbf,d3$logbf),
                     pr=list(d1$pr,d2$pr,d3$pr),
                     kappa=S,
                     p0=prior.bin.fn(0))
    d1$spp.models <- as.numeric(PPM$shared.pp[[1]])[-1]
    d1$mpp.models <- as.numeric(PPM$single.pp[[1]])[-1]
    d1$spp.str <- as.numeric(PPS$shared.pp[[1]])[-1]
    d1$mpp.str <- as.numeric(PPS$single.pp[[1]])[-1]
    return(list(d1=d1,d2=d2,d3=d3))
}


################################################################################

context("testing marginalpp")
## load data
testthat::skip_on_travis()
data <- read.table("/scratch/wallace/IL2RA/AD-AC-rep_100-bf-final.txt",header=TRUE,as.is=TRUE)

## analyse
data <- bothpp(data)
marg <- margpp(data)
adata <- expand.grid(1:nrow(marg$d1),1:nrow(marg$d2))
adata <- cbind(adata,marg$d1[adata$Var1,.(str,logbf,size)])
adata <- cbind(adata,marg$d2[adata$Var2,.(str,logbf,size)])
colnames(adata)[-c(1:2)] <- c("t1.str","t1.logbf","t1.size","t2.str","t2.logbf","t2.size")
adata$t1t2.logbf <- adata$t1.logbf + adata$t2.logbf - 46.6

## now perform same analysis as above
adata <- bothpp(adata)
amarg <- margpp(adata)
test_that("marginalpp2 works", {
    testthat::skip_on_cran()
    expect_equal(amarg$d1$mpp.orig,amarg$d1$mpp.direct)
    expect_equal(amarg$d1$mpp.orig,amarg$d1$mpp.new)
    expect_equal(amarg$d1$spp.orig,amarg$d1$spp.new)
})
## now all matches :)

## check same with marginalpp
test_that("marginalpp works for 2 diseases", {
    testthat::skip_on_cran()
    amarg2 <- margpp2(adata)
    expect_equal(amarg2$d1$spp.orig,amarg2$d1$spp.models)
    expect_equal(amarg2$d1$spp.orig,amarg2$d1$spp.str)
})

################################################################################
    ## make 3 way combinations
    adata <- expand.grid(1:nrow(marg$d1),1:nrow(marg$d2),1:nrow(marg$d2))
    adata <- cbind(adata,marg$d1[adata$Var1,.(str,logbf,size)])
    adata <- cbind(adata,marg$d2[adata$Var2,.(str,logbf,size)])
    adata <- cbind(adata,marg$d2[adata$Var3,.(str,logbf,size)])
    colnames(adata)[-c(1:3)] <- c("t1.str","t1.logbf","t1.size",
                                  "t2.str","t2.logbf","t2.size",
                                  "t3.str","t3.logbf","t3.size" )
    adata$t1t2.logbf <- adata$t1.logbf + adata$t2.logbf + adata$t3.logbf- 46.6
    ## now perform same analysis as above
    adata <- threepp(adata)
    amarg <- margthree(adata)

test_that("marginalpp works for 3 diseases", {
    testthat::skip_on_cran()
    expect_equal(amarg$d1$spp.orig,amarg$d1$spp.models)
    expect_equal(amarg$d1$spp.orig,amarg$d1$spp.str)
})


