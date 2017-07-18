##' Calculate marginal model posterior probabilities for each disease 
##'
##' @title Marginal PP for models sharing information between diseases
##' @param M1 model matrix for disease 1
##' @param M2 model matrix for disease 2
##' @param ABF1 log ABF for models in M1 for disease 1
##' @param ABF2 log ABF for models in M2 for disease 2
##' @param pr1 prior for models in M1
##' @param pr2 prior for models in M2
##' @param S single value or vector of values to consider for the
##'     sharing scale parameter
##' @param p0 prior probability of the null model
##' @return list of single.pp1 (pp for each model in M1 for disease 1,
##'     assuming ABF for all other models are approximately 1),
##' shared
##' @author Chris Wallace
marginalpp <- function(M, ABF, pr, kappa, p0) {
    #1, M2, M3, ABF1, ABF2, ABF3, pr1, pr2, pr3, S, p0) {
    ## are any of M1, M2 the null model?
    n <- length(M)
    if(n<2)
        stop("Need at least 2 diseases")
    if( length(ABF)!=n || length(pr)!=n )
        stop("M, ABF and pr need to have the same lengths")
    ncols <- sapply(M,function(x) dim(x)[2])
    if(sd(ncols)!=0)
        stop("all matrices in M must have the same SNPs (columns)")

    ## remove null model if included
    for(i in seq_along(M)) {
        wh <- which.null(M[[i]])
        if(length(wh)) {
            M[[i]] <- M[[i]][-wh,]
            ABF[[i]] <- ABF[[i]][-wh]
        }
    }

    ## unweighted pp
    pp <- mapply(function(pr1,ABF1) { calcpp(addnull(pr1,p0),addnull(ABF1,0)) },
                 pr, ABF, SIMPLIFY=FALSE)

    ## Q
    fun <- switch(n,
                  NULL,
                  "calcQ2",
                  "calcQ3")
    M <- lapply(M,t)
    if(is.null(fun))
        stop("calcQ not written for ",n," diseases yet")

    Q <- do.call(fun, c(M, lapply(pp,"[",-1)))

    ## alt prior
    app <- vector("list",n)
    for(i in seq_along(Q)) {
        tmp <- lapply(kappa, function(k) {
            a <- pr[[i]] * (1 + (k-1) * Q[[i]])
            a/sum(a)
        })
        tmp <- (1-p0) * do.call("cbind",tmp)
        M[[i]] <- addnull(M[[i]],0)
        app[[i]] <- calcpp(addnull(tmp,p0),addnull(ABF2,0))
    }

    maxpower <- n * (n-1) / 2
    app <- vector("list",n)
    for(i in seq_along(Q)) {
        tmp <- lapply(kappa, function(k) {
            s <- k^(1:maxpower)/maxpower
            a <- pr[[i]] * (1 + t((s-1) * t(Q[[i]])))
            rowSums(a)/sum(a)
        })
        tmp <- (1-p0) * do.call("cbind",tmp)
        M[[i]] <- addnull(M[[i]],0)
        app[[i]] <- calcpp(addnull(tmp,p0),addnull(ABF2,0))
    }

    M1 <- addnull(M1,0)
    app1 <- calcpp(addnull(altpr1,p0),addnull(ABF1,0)) 

    altpr2 <- lapply(S, function(s) {
        a <- pr2 * (1 + (s-1) * Q[[2]])
        a/sum(a)
    })
    M2 <- addnull(M2,0)
    altpr2 <- (1-p0) * do.call("cbind",altpr2) 
    app2 <- calcpp(addnull(altpr2,p0),addnull(ABF2,0)) 

    altpr3 <- lapply(S, function(s) {
        a <- pr3 * (1 + (s-1) * Q[[2]])
        a/sum(a)
    })
    M2 <- addnull(M2,0)
    altpr2 <- (1-p0) * do.call("cbind",altpr2) 
    app2 <- calcpp(addnull(altpr2,p0),addnull(ABF2,0)) 
    list(single.pp1=pp1,shared.pp1=app1,M1=M1,
         single.pp2=pp2,shared.pp2=app2,M2=M2,
         S=S)
}

which.null <- function(M) {
    rs <- rowSums(M)
    which(rs==0)
}

##' Calculate marginal model posterior probabilities for each disease 
##'
##' @title Marginal PP for models sharing information between diseases
##' @param M1 model matrix for disease 1
##' @param M2 model matrix for disease 2
##' @param ABF1 log ABF for models in M1 for disease 1
##' @param ABF2 log ABF for models in M2 for disease 2
##' @param pr1 prior for models in M1
##' @param pr2 prior for models in M2
##' @param S single value or vector of values to consider for the
##'     sharing scale parameter
##' @param p0 prior probability of the null model
##' @return list of single.pp1 (pp for each model in M1 for disease 1,
##'     assuming ABF for all other models are approximately 1),
##' shared
##' @export
##' @author Chris Wallace
marginalpp2 <- function(M1, M2, ABF1, ABF2, pr1, pr2, S, p0) {
    ## are any of M1, M2 the null model?
    rs1 <- rowSums(M1)
    if(any(rs1==0)) {
        wh <- which(rs1==0)
        M1 <- M1[-wh,]
        ABF1 <- ABF1[-wh]
    }
    rs2 <- rowSums(M2)
    if(any(rs2==0)) {
        wh <- which(rs2==0)
        M2 <- M2[-wh,]
        ABF2 <- ABF2[-wh]
    }
                 
    pp1 <- calcpp(addnull(pr1,p0),addnull(ABF1,0))
    pp2 <- calcpp(addnull(pr2,p0),addnull(ABF2,0))
    Q <- calcQ2(t(M1),t(M2),pp1[-1],pp2[-1])
    altpr1 <- lapply(S, function(s) {
        a <- pr1 * (1 + (s-1) * Q[[1]])
        a/sum(a)
    })
    altpr1 <- (1-p0) * do.call("cbind",altpr1)
    M1 <- addnull(M1,0)
    app1 <- calcpp(addnull(altpr1,p0),addnull(ABF1,0)) 

    altpr2 <- lapply(S, function(s) {
        a <- pr2 * (1 + (s-1) * Q[[2]])
        a/sum(a)
    })
    M2 <- addnull(M2,0)
    altpr2 <- (1-p0) * do.call("cbind",altpr2) 
    app2 <- calcpp(addnull(altpr2,p0),addnull(ABF2,0)) 
    list(single.pp1=pp1,shared.pp1=app1,M1=M1,
         single.pp2=pp2,shared.pp2=app2,M2=M2,
         S=S)
}
##' Internal function to add a null entry
##'
##' @title addnull
##' @param what vector or matrix
##' @param val value to add for the null model 
##' @return what with an additional first entry or first row filled with val
##' @author Chris Wallace
addnull <- function(what,val) {
    if(is.vector(what))
        return(c(val,what))
    rbind(matrix(val,1,ncol(what)),
          what)
}
               

##' Internal function to combine vector or matrix prior with vector of log ABF to get pp per model
##'
##' @title calcpp
##' @param pr vector or matrix of priors
##' @param lBF vector of log ABF
##' @return vector or matrix of pp
##' @author Chris Wallace
calcpp <- function(pr,lBF) {
    pp <- log(pr) + lBF
    if(is.matrix(pp)) {
        denom <- apply(pp,2,logsum)
        exp(t(pp) - denom)
    } else {
        return(exp(pp - logsum(pp)))
    }
}

##' Internal function, logsum (copied from coloc package)
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##'
##' ie, you want sum(x), but have x already stored in logs.  log(sum(exp(x))) might fail,
##' but logsum(x) should work.
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
##' @examples
##' x <- 1:10
##' log(sum(x))
##' MTFM:::logsum(log(x))
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

#' Internal function, logdiff
##'
##' This function calculates the log of the difference of the exponentiated
##' logs taking out the max, i.e. insuring that the difference is not negative
##'
##' ie you want log(exp(x) - exp(y))
##' @title logdiff
##' @param x numeric
##' @param y numeric
##' @return max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
##' @author Chris Wallace
##' @examples
##' x <- 1001:1010
##' y <- 1:10
##' log(x-y)
##' MTFM:::logdiff(log(x),log(y))
logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}

