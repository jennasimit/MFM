## library(annotSnpStats)
## library(snpStats)
## library(GUESSFM)
## library(mlogitBMA)
## library(BMA)
## library(Rcpp)
## library(RcppArmadillo)
## library(parallel)
## library(data.table)
## library(gtools)
## library(fields)
## library(dplyr)

## options(scipen=999)



#' @title Generate set of tag SNPs with LD r2 
#' @param r2 r-squared for tag SNP generation
#' @param G0 SnpMatrix object for control data
#' @param mydir file to save tags object to, if specified 
#' @export
make.tags.fn <- function(r2,G0,mydir=NA){
 tags <- tag(G0, tag.threshold = r2)
 if(!is.na(mydir)) {
  tfile <- file.path(mydir,paste("tags-r2-",r2,".RData",sep="") )
  message("saving tags object to ", tfile)
  save(tags, file = tfile)
  }
 return(tags)
}

####


#' @title Find "best" tag SNPs for a trait to carry forward in models
#' @param snpG is a "SnpMatrix" object
#' @param y is the phenotype vector
#' @param tags is a tags object
#' @param mppthr is the threshold for selecting the "best" SNPs with marginal PP > mppthr
#' @param mydir is a directory name where to save stochastic search results on the tag SNPs
#' @return data.frame of tag SNPs with MPP > mppthr after running GUESSFM on tag SNPs
#' @export
gfm.sel.snps.fn <- function(snpG,y,tags,mppthr=0.001,mydir) {
 snpG
 tsnp <- snpG[,unique(tags@tags)]  
 run.bvs(X=tsnp,Y=y,nexp=3,tag.r2=NA,nsave=1000,gdir=mydir,wait=TRUE) 
 d <- read.snpmod(mydir)
 bestsnps <- best.snps(d,pp.thr=mppthr)
return(bestsnps)
}


format.mod.fn <- function(k,out) {
#' called by T1mods.fn
 ind <- which(out[k,]==1)
 if(length(ind)>0) {
  mod <- paste(names(out[k,][ind]),sep="",collapse="%")
  } else {mod <- ""}
  return(mod)
	}


#' @title Generate all possible models for model sizes in a pre-specified range
#' @param mT1 is the minimum model size
#' @param MT1 is the maximum model size
#' @param msnps is the vector of SNPs that will be selected for models
#' @return data.frame of models (e.g. snp1%snp2 for a 2-SNP model of snp1 and snp2) and their sizes 
#' @export
#' @author Jenn Asimit
T1mods.fn <- function(mT1,MT1,msnps) {
 s <- length(msnps)
 if(mT1>0) mc <- mT1:MT1
 if(mT1==0) mc <- (mT1+1):MT1 # start at 1 instead of 0 
 modT1 <- vector("list",length(mc))
  for(i in 1:length(mc)) {
   mci <- combn(msnps,mc[i],simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s,dimnames=list(NULL,msnps))
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],msnps)] <- 1 
  }
 #' if(mT1 > 0) T1mods <- NULL
 #' if(mT1 == 0) T1mods <- matrix(rep(0,s),nrow=1)
 T1mods <- NULL
 for(i in 1:length(mc)) T1mods <- rbind(T1mods,modT1[[i]])
 
 j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
 t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
 t1size <- apply(T1mods,1,sum)
 T1mod <- data.frame(mod=t1mod,size=t1size,stringsAsFactors=FALSE)
 
 return(T1mod)
}

####

check.fn <- function(k,msep,out,gnames) {
#' called by MPP.fn
     g <- length(gnames)
     p1 <- numeric(g) 
     for(j in 1:g) { 
     ind1 <- gnames[j] %in% msep[[k]]
     if(ind1) p1[j] <- out[k] 
     }
     return(p1)
    	}

####

sep.fn <- function(k,mnames) {
#' called by MPP.fn
   msep <- unlist(strsplit(as.character(mnames[k]),"%"))
   return(msep)
   }

####

#' @title Find marginal PP of inclusion for SNPs (MPPi) in marginal models at each kappa value
#' @param PP1 is a matrix of marginal model PP for a trait (i.e. one list component of output from sharedmPP.fn),
#' @return a matrix of MPPi, where rows are SNPs and columns are kappa
#' @export
#' @author Jenn Asimit
MPP.fn<-function(PP1) {
 mnames <- rownames(PP1) 
  msep <- apply(matrix(1:length(mnames),ncol=1),1,sep.fn,mnames)
  gnames <- unique(unlist(msep)) # snps     
   mpp1 <- NULL
   for(k in 1:dim(PP1)[2]) {
    tmp1 <- apply(matrix(1:length(mnames),ncol=1),1,check.fn,msep,PP1[,k],gnames)  
    mpp1 <- rbind(mpp1,apply(tmp1,1,sum) )
    }    
   mpp1 <- data.frame(mpp1,row.names=colnames(PP1))  
   names(mpp1)<-gnames    
return(t(mpp1))
}

####

#' @title Expand tag SNP models
#' @param best is the output from T1mods.fn, a data.frame of the models and their sizes
#' @param tags is the set of tag SNPs; a tags object
#' @return a list with three components mod (tag SNP model), size (model size), str (expanded tag SNP model)
#' @export
expand.tags.bf <- function(best, tags) {
   
    bsnps <- unique(unlist(strsplit(as.character(best$mod),"%")))    
    check <- which(bsnps=="0" | bsnps=="1") 
    if(length(check)>0) bsnps <- bsnps[-check]
    B <- dim(best)[1]
    wh <- which(make.names(tags(tags)) %in% bsnps)
           
    
    if (!length(wh)) 
        stop("none of the supplied tags are amongst the best SNPs in d")
    proxies <- split(make.names(snps(tags)[wh]), make.names(tags(tags)[wh]))
   
    if (!all(bsnps %in% names(proxies))) 
        stop("not all model SNPs found in tags object")
    message("expanding tags for ", B, " models over ", length(proxies), 
        " tag SNPs, tagging a total of ", length(unlist(proxies)), 
        " SNPs.")
    pm <- mclapply(as.list(1:B), function(i) {
        if (best[i, "size"] == 0) {
            pm.str <- ""
        }
        else {
            pm.str <- apply(do.call(expand.grid, proxies[unlist(strsplit(as.character(best$mod[i]),"%"))]), 
                1, makestr)
        }
        gc()
        
        return(pm.str)
    })
    npm <- sapply(pm, length)
    neighb <- as.data.table(best)[rep(1:length(pm), times = npm), 
        ]
    setnames(neighb, sub("str", "index.str", names(neighb)))
    neighb[, `:=`(str, unlist(pm))]
    best <- as.data.frame(neighb,stringsAsFactors=FALSE)
   return(best)
}

####

#' @title Find best tag SNP models, then expand and calculate joint PP at these expanded models
#' @param kappa is a vector of sharing scale values for the joint prior
#' @param snp.data.list is a list of SnpMatrix objects for each trait
#' @param pheno.list is a list of phenotype vectors for each trait
#' @param tags tags object of SNPs from common controls 
#' @param mppthr threshold for "best" SNPs marginal PP of inclusion
#' @param traits is a vector of trait names
#' @param mydir is a directory stub to save marginal stochastic search (GUESSFM) results on tag SNPs 
#' @return a list with components pp1, pp2 (each is a data.frame of the marginal model PPs at each kappa value)
#' @export
#' @author Jenn Asimit
sharedmPP.fn <- function(kappa,snp.data.list,pheno.list,tags,mppthr,traits,mydir) {
K <- length(traits) 
mydirs <- paste(mydir,1:K,sep="_")

#' find "best" tag snps
msnps <- c()
for(k in 1:K) {
 tsnps <- gfm.sel.snps.fn(snpG=snp.data.list[[k]],y=pheno.list[[k]],tags=tags,mppthr=mppthr,mydir=mydirs[k])
 msnps <- union(msnps,tsnps[,"var"])
 }
 

#' find min and max model sizes for each trait
mT <- numeric(K)
MT <- numeric(K)
for(k in 1:K) {
d <- read.snpmod(mydirs[k])
dx <- expand.tags(d,tags)
 best <- best.models(dx,pp.thr=0.0001)
 abf <- abf.calc(y=pheno.list[[k]],x=snp.data.list[[k]] ,models=best$str,family="binomial") 
sm <- abf2snpmod(abf,expected=3,nsnps=605)
modsizes <- pp.nsnp(sm,expected=3)
srange <- as.numeric(names(modsizes$trait)[which(modsizes$trait>0.001)])
mT[k] <- min(srange)
MT[k] <- max(srange) 
}

   
#' generate all models from "best" tag snps
Tmods <- vector("list",K)
for(k in 1:K) Tmods[[k]] <- T1mods.fn(mT[k],MT[k],msnps) 
 
 
 G.mat <- vector("list",length(traits))
 DATA <- vector("list",length(traits))
 for(k in 1:K) {
 G.mat[[k]] <- as(snp.data.list[[k]],"numeric")
 DATA[[k]] <- data.frame(Y=pheno.list[[k]],G.mat[[k]])
}
 
    
 e.list <- vector("list",K)
 for(k in 1:K) e.list[[k]] <- expand.tags.bf(Tmods[[k]],tags)


 ABF <- vector("list",K)
 for(k in 1:K) ABF[[k]] <- abf.calc(y=DATA[[k]][,1],x=DATA[[k]][,-1],models=e.list[[k]]$str,family="binomial")[[1]]$lBF 

 
 STR <- vector("list",K)
 for(k in 1:K) STR[[k]] <- e.list[[k]]$str
 pr <- vector("list",K)
 for(k in 1:K) pr[[k]] <- prior.bin.fn(e.list[[k]]$size)
 
 p0 <- prior.bin.fn(0)


 mpp <- marginalpp(STR, ABF, pr, kappa, p0)

 pp <- vector("list",length(traits))
 for(k in 1:length(traits)) {
  pp[[k]] <- t(mpp$shared.pp[[k]])
 rownames(pp[[k]]) <- mpp$STR[[k]]
 colnames(pp[[k]]) <- paste("pp",kappa,sep=".") 
 }
   
 return(pp)
}

