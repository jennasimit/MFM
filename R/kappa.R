#' @title Prior probability (in log units) for single trait models
#' @param x Vector of number of SNPs in each model
#' @param n Number of SNPs in region
#' @param v Expected number of causal variants for region
#' @return Vector of log prior probabilities for each single trait model of given size 
p.fn <- function(x,n,v) dbinom(x,size=n,prob=v/n,log=TRUE)


#' @title Kappa selection
#' @param Nd Number of diseases to fine-map
#' @param Ps2 Probability of any shared causal variant between a disease pair
#' @param n Number of SNPs in region
#' @return Optimal kappa value
#' @export
kappaOpt <-function(Nd,Ps2,nsnps) {
 Ps <- 1-(1-Ps2)^(Nd-1) # Pr(given disease shares a variant with at least one other disease)
 pr <- p.fn(0:nsnps,nsnps,2)
 if(Nd==2) kappa.opt <- kappa2(Ps,nsnps,100,pr)
 if(Nd==3) kappa.opt <- kappa3(Ps,nsnps,100,pr)
 if(Nd==4) kappa.opt <- kappa4(Ps,nsnps,10,pr)
 return(kappa.opt)
}


os <- function(kappa,prob) {
    n <- length(prob)-1
    pij <- function(i,j) {
        ## choose(n-i,j) * prob[i] * prob[j] / ( (choose(n,j) - choose(n-i,j))*kappa + choose(n-i,j) )
        denom <- 1 + kappa * (choose(n,j) / choose(n-i,j) - 1)
        p <- prob[i+1] * prob[j+1] / denom
    }
    pn <- unlist(lapply(0:n, function(i) {
        unlist(lapply(0:n, function(j) pij(i,j)))
    }))
    pn <- sum(pn)
    log(1-pn) - log(pn)
}
##' Calculate kappa for a given region and priors on the number of causal variants in the region, and the prior odds of causal variant sharing
##'
##' @title calckappa
##' @return the value of kappa corresponding to supplied target.odds
##' @author Chris Wallace
##' @export
##' @param nsnps number of snps in the region
##' @param p probability any snp is causal (assumes a binomial distribution for the number of causal snps, with expected number = p * nsnps)
##' @param target.odds prior odds of any sharing of causal variants between any pair of traits
calckappa <- function(nsnps,p,target.odds) {
    prob <- dbinom(0:nsnps,size=nsnps,prob=p)
    f <- function(kappa) {
        abs(odds_sharing(kappa,prob) - log(target.odds))
    }
    ## k <- seq(1,10,by=0.1)
    ## o <- sapply(k,odds_sharing,p=prob)
    ## plot(k,o,type="b")
    optimize(f, c(1,1000))$minimum
}
