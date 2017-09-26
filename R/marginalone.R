##' Calculate marginal model posterior probabilities for one disease
##' incorporating information from others
##'
##' Given a list of model matrices and log ABFs, this function
##' calculates the marginal model posterior probabilities for thr
##' FIRST disease without ever calculating the joint Bayes Factors for
##' all cross-disease model configurations, which would require large
##' amounts of memory.
##'
##' This aims to be faster than marginalpp by ignoring contributions
##' from other diseases to models with individual PP < eps.  This is
##' an approximation, and the effects of having too large an eps may
##' induce inaccuracies.  However, if you want to pull in information
##' from several diseases, this may be the only way to do it within
##' achievable computer time.
##'
##' @title Marginal PP for models sharing information between diseases
##' @param STR list of models for diseases 1, 2, ..., n, each given in
##'     the form of a character vector, with entries
##'     \code{"snp1\%snp2\%snp3"}. The null model is given by
##'     \code{"1"} OR \code{"0"}.  It is assumed that all elements of
##'     ABF, PP and pr below follow this same order.
##' @param ABF list of log(ABF) vectors for diseases 1, 2, ...
##' @param PP list of posterior probability vectors for diseases 1, 2,
##'     ...
##' @param pr list of prior probabilities for the models in M
##' @param kappa single value or vector of values to consider for the
##'     sharing scale parameter.  the value of kappa=1 must be
##'     included, and if not will be prepended.
##' @param p0 prior probability of the null model
##' @param eps models for all but the first disease are dropped if
##'     their PP < eps.  Ie set PP=0 if PP<eps.  This is an
##'     APPROXIMATION, eps should be as small as your computing
##'     facilities allow.
##' @return list of: - single.pp: list of pp for each model in
##'     STR[[i]] for disease i - shared.pp: list of pp for each model
##'     in STR[[i]] for disease i, - STR: not quite as input,
##'     reordered so null model is first row - ABF: not quite as
##'     input, repordered so null model is first row - kappa: as
##'     supplied
##' @export
##' @author Chris Wallace
marginalone <- function(STR, ABF, PP, pr, kappa, p0, tol=0.0001, eps=1e-6) {
    n <- length(STR)
    if(n<2)
        stop("Need at least 2 diseases")
    if( length(ABF)!=n || length(pr)!=n | length(PP)!=n )
        stop("STR, ABF, PP and pr need to have the same lengths")
    SS <- lapply(STR,strsplit,"%")
    SS <- lapply(SS,setdiff,c("0","1"))
    usnps <- sort(unique(unlist(SS)))
    if(!(1 %in% kappa))
        kappa <- c(1,kappa)
    
    ## remove null model if included
    PP.nonull <- PP
    for(i in seq_along(STR)) {
        wh <- which(STR[[i]] %in% c("0","1"))
        if(length(wh)) {
            STR[[i]] <- STR[[i]][-wh]
            ABF[[i]] <- ABF[[i]][-wh]
            pr[[i]] <- pr[[i]][-wh]
            PP.nonull[[i]] <- PP[[i]][-wh]
            ## PP[[i]] <- c(PP[[i]][wh],PP[[i]][-wh])
            PP[[i]] <- PP[[i]][-wh]
            ## ABF[[i]] <- addnull(ABF[[i]],0)
            ## pr[[i]] <- addnull(pr[[i]],p0)
            ## PP[[i]] <- addnull(PP[[i]], ABF[[i]][1] * pr[[i]][1] / sum(ABF[[i]] * pr[[i]]))
        } 
        PP[[i]] <- addnull(PP[[i]], calcpp(addnull(pr[[i]],p0), addnull(ABF[[i]],0))[1])
    }
    
    STR.i <- lapply(SS, function(ss) {
        lapply(ss,function(x) as.integer(factor(x,levels=usnps)))
    })
    names(STR.i) <- NULL
    
    ## unweighted pp - as input
    ## pp <- mapply(function(pr1,ABF1) {
    ##     calcpp(addnull(pr1,p0),addnull(ABF1,0)) },
    ##     pr, ABF, SIMPLIFY=FALSE)
    names(PP.nonull) <- NULL

    ## trim subsequent diseases
    message("trimming PP < ",eps," for diseases 2..",n)
    message("initial lengths: ",paste(lapply(PP,length,collapse=", ")))
    PP[-1] <- lapply(PP[-1], function(x) x[ which(x) >= eps ])
    message("trimmed lengths: ",paste(lapply(PP,length,collapse=", ")))
    
    ## Q
    fun <- switch(n,
                  NULL,
                  "calcQone2",
                  "calcQone3",
                  "calcQone4")
    if(is.null(fun))
        stop("calcQone not written for ",n," diseases yet")

    Q <- do.call(fun, c(STR.i, PP.nonull)) #lapply(pp,"[",-1)))
    
    ## alt prior
    maxpower <- n * (n-1) / 2
    tmp <- lapply(kappa, function(k) {
        if(n==2) {
            a <- pr[[1]] * (1 + (k-1) * Q)
        } else {
            s <- k^((1:maxpower)/maxpower)
            a <- pr[[1]] * (1 + colSums((s-1) * t(Q)))
        }
        a#/sum(a)
    })
    tmp <- do.call("cbind",tmp)
    alt.prior <- addnull(tmp,p0)
    alt.pp <- calcpp(alt.prior,addnull(ABF[[1]],0))
    pr <- lapply(pr, addnull, p0)
    STR[[1]] <- addnull(STR, "1")
    alt.pp <- t(alt.pp)
    
    ## checks
    wh <- which(kappa==1)
    sumsq <- sum((PP[[1]] - alt.pp)^2)
    if((sumsq>tol)) {
        warning("trait ",1," kappa=1 PP does not match input PP, sumsq=",sumsq,
                "which is > tol.\nsuggests you need to include more models in the calculation")
    }
    
    list(single.prior=pr[[1]], single.pp=PP[[1]],
         shared.prior=alt.prior,shared.pp=alt.pp,
         STR=STR[[1]],kappa=kappa)  
}
    

