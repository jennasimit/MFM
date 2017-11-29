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
 Ps <- 1-Ps2^(Nd-1) # Pr(given disease shares a variant with at least one other disease)
 pr <- p.fn(0:nsnps,nsnps,2)
 if(Nd==2) kappa.opt <- kappa2(Ps,nsnps,100,pr)
 if(Nd==3) kappa.opt <- kappa3(Ps,nsnps,100,pr)
 if(Nd==4) kappa.opt <- kappa4(Ps,nsnps,10,pr)
 return(kappa.opt)
}


