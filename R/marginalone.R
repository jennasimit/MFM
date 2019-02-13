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
##' @param pr list of prior probabilities for the models in STR
##' @param kappa single value or vector of values to consider for the
##'     sharing scale parameter.  the value of kappa=1 must be
##'     included, and if not will be prepended.
##' @param p0 prior probability of the null model
##' @param fthr models for all but the first disease are retained if
##'     their cumsum(PP) < fthr.  Ie set \code{PP[j]}=0 if SNP j is not in the smallest set of SNPs that satisfy cumsum(PP) < fthr.  This is an
##'     APPROXIMATION, eps should be as close to 1 as your computing
##'     facilities allow.
#' @param N0 number of shared controls
#' @param ND list of number of cases for a set of diseases
##' @return list of: - single.pp: list of pp for each model in
##'     \code{STR[[i]]} for disease i - shared.pp: list of pp for each model
##'     in \code{STR[[i]]} for disease i, - STR: not quite as input,
##'     reordered so null model is first row - ABF: not quite as
##'     input, repordered so null model is first row - kappa: as
##'     supplied
##' @export
##' @author Chris Wallace
marginalone <- function(STR, ABF, PP, pr, kappa, p0, fthr=0.99,N0,ND) {
    n <- length(STR)
    if (n < 2) 
        stop("Need at least 2 diseases")
    if (length(ABF) != n || length(pr) != n | length(PP) != n) 
        stop("STR, ABF, PP and pr need to have the same lengths")

dis <- names(STR)
N <- sum(unlist(ND))+N0
Mk <- vector("list",n)
 PP.adj <- vector("list",n) # adjusted PP for computations
 pr.adj <- vector("list",n)
 for(j in 1:n) {
  Mk[[j]] <- unlist(lapply(strsplit(STR[[j]],"%"),length)) # model sizes
  eta <- exp(Mk[[j]]*.5*log((ND[[j]]+N0)/N))
  PP.adj[[j]] <- PP[[j]]*eta/sum(PP[[j]]*eta)
  pr.adj[[j]] <- pr[[j]]*eta/sum(pr[[j]]*eta)
   }
  names(PP.adj) <- dis
  names(pr.adj) <- dis

	print(paste("trimming cumsumPP > ", fthr, " for diseases 2..", n))
    print(paste("initial lengths: ", paste(sapply(PP, length), collapse = ", ")))
   
   PP0 <- PP
   pr0 <- pr
   pr[[1]] <- pr.adj[[1]]
   PP[[1]] <- PP.adj[[1]]
   
	for (i in 2:n) {
	    ind <- order(PP[[i]],decreasing=TRUE)
	    trun <- min(which(cumsum(sort(PP[[i]],decreasing=TRUE))>fthr))	    
        wh <- ind[1:trun]
        if (length(wh)) {
            STR[[i]] <- STR[[i]][wh]
            ABF[[i]] <- ABF[[i]][wh]
            pr[[i]] <- pr.adj[[i]][wh]
            PP[[i]] <- PP.adj[[i]][wh]
            pr0[[i]] <- pr0[[i]][wh]
            PP0[[i]] <- PP0[[i]][wh]
        } else {
            pr[[i]] <- pr.adj[[i]]
            PP[[i]] <- PP.adj[[i]]
 				}       
    }
   
    print(paste("trimmed lengths: ", paste(sapply(PP, length), collapse = ", ")))

    
    SS <- lapply(STR, strsplit, "%")
    SS <- lapply(SS, setdiff, c("0", "1"))
    usnps <- sort(unique(unlist(SS)))
    if (!(1 %in% kappa)) 
        kappa <- c(1, kappa)
    PP.nonull <- PP
    for (i in seq_along(STR)) {
        wh <- which(STR[[i]] %in% c("0", "1"))
        if (length(wh)) {
            STR[[i]] <- STR[[i]][-wh]
            ABF[[i]] <- ABF[[i]][-wh]
            pr[[i]] <- pr[[i]][-wh]
            pr0[[i]] <- pr0[[i]][-wh]
            PP.nonull[[i]] <- PP[[i]][-wh]
            PP[[i]] <- PP[[i]][-wh]
            PP0[[i]] <- PP0[[i]][-wh]
        }
        PP[[i]] <- addnull(PP[[i]], calcpp(addnull(pr[[i]], p0), 
            addnull(ABF[[i]], 0))[1])
        PP0[[i]] <- addnull(PP0[[i]], calcpp(addnull(pr0[[i]], p0), 
            addnull(ABF[[i]], 0))[1])    
    }
      STR.i <- lapply(SS, function(ss) {
        lapply(ss, function(x) as.integer(factor(x, levels = usnps)))
    })
    names(STR.i) <- NULL
    names(PP.nonull) <- NULL
  
    
    fun <- switch(n, NULL, "calcQone2", "calcQone3", "calcQone4")
    if (is.null(fun)) 
        stop("calcQone not written for ", n, " diseases yet")
    
      
    Q <- do.call(fun, c(STR.i, PP.nonull))
    maxpower <- n * (n - 1)/2
    tmp <- lapply(kappa, function(k) {
        if (n == 2) {
            a <- pr[[1]] * (1 + (k - 1) * Q)
        }
        else {
            s <- k^((1:maxpower)) #/maxpower)
            a <- pr[[1]] * (1 + colSums((s - 1) * t(Q)))
        }
        a
    })
    tmp <- do.call("cbind", tmp)
    alt.prior <- addnull(tmp, p0)
    alt.pp <- calcpp(alt.prior, addnull(ABF[[1]], 0))
    pr <- lapply(pr, addnull, p0)
    STR[[1]] <- addnull(STR[[1]], "1")
    alt.pp <- t(alt.pp)
   # wh <- which(kappa == 1)
   # sumsq <- sum((PP0[[1]] - alt.pp[,wh])^2) 
   # if ((sumsq > tol)) {
   #     warning("trait ", 1, " kappa=1 PP does not match input PP, sumsq=", 
    #        sumsq, "which is > tol.\nsuggests you need to include more models in the calculation")
    #}
    list(single.prior = pr[[1]], single.pp = PP[[1]], shared.prior = alt.prior, 
        shared.pp = alt.pp, STR = STR[[1]], kappa = kappa)
}

