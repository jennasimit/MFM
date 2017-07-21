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



#' Generate set of tag SNPs from control genotype matrix G0 ("SnpMatrix" object for controls set)
#' @title Generate set of tag SNPs with LD r2 
#' @param r2 r-squared for tag SNP generation
#' param G0 SnpMatrix object for control data
#' param mydir file to save tags object to, if specified 
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

#' find tag SNPs with MPP > mppthr after running GUESSFM on tag SNPs and save to file
#' print table of posterior distribution for number of causal variants in GFM models
#' @export
gfm.sel.snps.fn <- function(snpG,y,tags,mppthr=0.001,mydir) {
#' snpG is a "SnpMatrix" object; y is the phenotype vector; mppthr is the threshold for selecting
#' the "best" SNPs with marginal PP > mppthr
 snpG
 tsnp <- snpG[,unique(tags@tags)]  
 run.bvs(X=tsnp,Y=y,nexp=3,tag.r2=NA,nsave=1000,gdir=mydir,wait=TRUE) 
 d <- read.snpmod(mydir)
 bestsnps <- best.snps(d,pp.thr=mppthr)

 print(pp.nsnp(d,expected=3))

return(bestsnps)
}

####

format.mod.fn <- function(k,out) {
#' called by T1mods.fn
 ind <- which(out[k,]==1)
 if(length(ind)>0) {
  mod <- paste(names(out[k,][ind]),sep="",collapse="%")
  } else {mod <- ""}
  return(mod)
	}

####

T1mods.fn <- function(mT1,MT1,msnps) {
 #' mT1 and MT1 are the min and max number of causals in a model for trait t1
 #' msnps is the vector of snps to consider in the models
 s <- length(msnps)
 if(mT1>0) mc <- mT1:MT1
 if(mT1==0) mc <- (mT1+1):MT1 # start at 1 instead of 0 and do 0 snp model later
 modT1 <- vector("list",length(mc))
  for(i in 1:length(mc)) {
   mci <- combn(msnps,mc[i],simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s,dimnames=list(NULL,msnps))
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],msnps)] <- 1 
  }
 if(mT1 > 0) T1mods <- NULL
 if(mT1 == 0) T1mods <- matrix(rep(0,s),nrow=1)
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

#' @title Calculate joint PP at tag SNP models and then expanded models
#' @param t1snp.data a SnpMatrix object for trait 1
#' @param t1pheno a phenotype vector for trait 1
#' @param t2snp.data a SnpMatrix object for trait 2
#' @param t2pheno a phenotype vector for trait 2
#' @param tags tag SNPs from common controls 
#' @param mppthr threshold for "best" SNPs MPP
#' @param mT1 minimum model size for trait 1
#' @param MT1 maximum model size for trait 1
#' @param mT2 minimum model size for trait 2
#' @param MT2 maximum model size for trait 2
#' @return a list with components pp1, pp2 (each is a data.frame of the marginal model PPs at each kappa value)
#' @export
sharedmPP.fn <- function(kappa,t1snp.data,t1pheno,t2snp.data,t2pheno,
			mppthr,mT1,MT1,mT2,MT2,
			trait1="T1",trait2="T2") {

#' find "best" tag snps
  t1snps <- gfm.sel.snps.fn(snpG=t1snp.data,y=t1pheno,tags=tags,mppthr=mppthr,mydir)
  t2snps <- gfm.sel.snps.fn(snpG=t2snp.data,y=t2pheno,tags=tags,mppthr=mppthr,mydir)
  msnps <- union(t1snps[,"var"],t2snps[,"var"])
     s <- length(msnps)
#' NEED to automate MT1,MT2 selection
   
#' generate all models from "best" tag snps
 T1mod <- T1mods.fn(mT1,MT1,msnps)
 T2mod <- T1mods.fn(mT2,MT2,msnps)
 
 G1.mat <- as(t1snp.data,"numeric")
 data1 <- data.frame(Y=t1pheno,G1.mat)

 G2.mat <- as(t2snp.data,"numeric")
 data2 <- data.frame(Y=t2pheno,G2.mat)

 
    D1 <- data.frame(mod=T1mod$mod,size=T1mod$size)
    D2 <- data.frame(mod=T2mod$mod,size=T2mod$size)

 e1 <- expand.tags.bf(D1,tags)
 e2 <- expand.tags.bf(D2,tags)
 
 bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=e1$str,family="binomial")[[1]] 
 bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=e2$str,family="binomial")[[1]] 
 
 
 STR <- list(e1$str,e2$str)
 ABF <- list(bf1$lBF,bf2$lBF)
 pr1 <- prior.bin.fn(e1$size)
 pr2 <- prior.bin.fn(e2$size)
 pr <- list(pr1,pr2)
 p0 <- prior.bin.fn(0)

 mpp <- marginalpp(STR, ABF, pr, kappa, p0)

 pp1 <- t(mpp$shared.pp[[1]])
 rownames(pp1) <- mpp$STR[[1]]
 colnames(pp1) <- paste("pp",kappa,sep=".")
  
 pp2 <- t(mpp$shared.pp[[2]])
 rownames(pp2) <- mpp$STR[[2]]
 colnames(pp2) <- paste("pp",kappa,sep=".")
  
 return(list(pp1=pp1,pp2=pp2))
}
  	 
#### these should be deleted

#' @title Expand tag models and calculate joint BFs at expanded models
#' @param ppbf output from PPBF.tags.fn
#' @param tags set of tag SNPs from common controla
#' @param t1snp.data a SnpMatrix object for trait 1
#' @param t1pheno a phenotype vector for trait 1
#' @param t2snp.data a SnpMatrix object for trait 2
#' @param t2pheno a phenotype vector for trait 2
#' @param trait1 name of trait 1, for filename purposes
#' @param trait2 name of trait 2, for filename purposes
#' @export
BFexpanded.fn <- function(ppbf,tags,t1snp.data,t1pheno,t2snp.data,t2pheno,
			mydir,trait1,trait2) {
 pp <- ppbf$pp
 logBF0 <- ppbf$logBF0
 
 indbest <- which(pp$PP>0.00001)
 BFkeep <- pp[-indbest,] # expand low pp models to see how many exist for each model and each of these will have the same bf as the tag model

 BFexp.notbest <- T1T2bfexp.low.fn(BFkeep,tags,logBF0)

#' expand tags at all marginal models that are part of a "best" joint model
 tmp <- data.frame(logbf=pp$t1.logbf,mod=pp$t1.mod,size=pp$t1.size)
 best1 <- unique(tmp[indbest,])

 tmp <- data.frame(logbf=pp$t2.logbf,mod=pp$t2.mod,size=pp$t2.size)
 best2 <- unique(tmp[indbest,])

e1 <- expand.tags.bf(best1,tags)
e2 <- expand.tags.bf(best2,tags)

G1.mat <- as(t1snp.data,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)

G2.mat <- as(t2snp.data,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)

bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=e1$str,family="binomial") 
bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=e2$str,family="binomial") 

T1bfexp <- e1
T1bfexp$logbf <- bf1[[1]][,"lBF"]

T2bfexp <- e2
T2bfexp$logbf <- bf2[[1]][,"lBF"]

#' at expanded tags, re-fit joint models 
T1T2bfexp <- T1T2bfexp.best.fn(T1T2abf,T1bfexp,T2bfexp,indbest,logBF0)

#' merge expanded re-fit joint with expanded non-re-fit joint 
tmp<-BFexp.notbest[,c("t1.logbf","t1.str","t1.size","t2.logbf","t2.str","t2.size","t1t2.logbf","t1t2.Nexpmod")]
BF <- rbind(T1T2bfexp,tmp)

traits <- paste(trait1,"-",trait2,sep="")
fname <- file.path(mydir,paste(traits,"-bf-final.txt",sep=""))
write.table(BF,fname,quote=FALSE,row.names=FALSE,col.names=TRUE)

return(BF) 
}

###

#' @title Calculate joint PP at a set of sharing scales
#' @param shared vector of sharing scales 
#' @param bf output from BFexpanded.fn
#' @param tags set of tag SNPs for control set
#' @param trait1 name of trait 1, for filename purposes
#' @param trait2 name of trait 2, for filename purposes
#' @mydir directory to save output to
#' @return a list of PPmarg (PP for each trait) and mpp (marginal PP for each trait) evalauted at each sharing scale value
#' @export
PPshared.fn <- function(shared,bf,tags,trait1,trait2,mydir) {

 ns <- length(shared)
 PP <- NULL
 for(k in 1:ns) PP <- rbind(PP,PP.fn(bf,shared=shared[k],s=length(snps(tags)),mT1=3,mT2=3,details=FALSE))
 row.names(PP) <- shared
 traits <- paste(trait1,"-",trait2,sep="")
 fname <- file.path(mydir,paste(traits,"-PP-shared.txt",sep=""))
 PP <- t(PP)
 write.table(PP,fname,row.names=TRUE,col.names=TRUE,quote=FALSE)

 PPall <- data.frame(PP,BF$t1t2.Nexpmod)

 alltraits <- c(trait1,trait2)
 K <- length(alltraits)
 
 mpp <- vector("list",K)
 PPmarg <- PP.marg.fn(PPall,K)
 for(k in 1:K) {
 fname <- file.path(fdir,paste(alltraits[k],"-PP-shared-",traits,".txt",sep=""))
 write.table(PPmarg[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
 mpp[[k]] <- MPP.fn(PPmarg[[k]])
 fname <- file.path(fdir,paste(alltraits[k],"-MPP-shared-",traits,".txt",sep=""))
 write.table(mpp[[k]],fname,row.names=TRUE,col.names=TRUE,quote=FALSE)
  }
  return(list(PPmarg=PPmarg,mpp=mpp))
}

####

#' @title Calculate joint BFs and joint PP at tag SNP models
#' @param t1snp.data a SnpMatrix object for trait 1
#' @param t1pheno a phenotype vector for trait 1
#' @param t2snp.data a SnpMatrix object for trait 2
#' @param t2pheno a phenotype vector for trait 2
#' @param tags tag SNPs from common controls 
#' @param mppthr threshold for "best" SNPs MPP
#' @param mT1 minimum model size for trait 1
#' @param MT1 maximum model size for trait 1
#' @param mT2 minimum model size for trait 2
#' @param MT2 maximum model size for trait 2
#' @return a list with components pp (data.frame of joint PP, prior, logprior, logBF, and for each trait, model, logBF, model size) and logBF0, the offset term needed to calculate joint BFs
#' @export
PPBF.tags.fn <- function(t1snp.data,t1pheno,t2snp.data,t2pheno,
			tags,mppthr,mT1,MT1,mT2,MT2,
			trait1="T1",trait2="T2") {

#' find "best" tag snps
 t1snps <- gfm.sel.snps.fn(snpG=t1snp.data,y=t1pheno,tags=tags,mppthr=mppthr)
 t2snps <- gfm.sel.snps.fn(snpG=t2snp.data,y=t2pheno,tags=tags,mppthr=mppthr)
 msnps <- union(t1snps[,"var"],t2snps[,"var"])
 s <- length(msnps)

#' generate all models from tag snps
 T1mod <- T1mods.fn(mT1,MT1,msnps)
 T2mod <- T1mods.fn(mT2,MT2,msnps)

 G1.mat <- as(t1snp.data,"numeric")
 data1 <- data.frame(Y=t1pheno,G1.mat)

 G2.mat <- as(t2snp.data,"numeric")
 data2 <- data.frame(Y=t2pheno,G2.mat)

 bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 
 bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 

 bf1$size <- t1size
 bf2$size <- t2size

 logBF0 <- logBF0.fn(msnps,t1pheno,G1.mat,t2pheno,G2.mat,trait1,trait2)

 T1T2abf <- T1T2abfexp.fn(bf1,bf2,logBF0)

 pp<-PP.fn(T1T2abf,shared=1,s=length(unique(tags(tags))),mT1=3,mT2=3)
 
 return(list(pp=pp,logBF0=logBF0))
}
