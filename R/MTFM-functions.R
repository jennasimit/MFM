library(annotSnpStats)
library(snpStats)
library(GUESSFM)
library(mlogitBMA)
library(BMA)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(data.table)
library(gtools)
library(fields)
library(dplyr)

options(scipen=999)



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
gfm.sel.snps.fn <- function(snpG,y,tags,mppthr=0.001) {
#' snpG is a "SnpMatrix" object; y is the phenotype vector; mppthr is the threshold for selecting
#' the "best" SNPs with marginal PP > mppthr
 snpG
 tsnp <- snpG[,unique(tags@tags)]  
 run.bvs(X=tsnp,Y=y,nexp=3,tag.r2=NA,nsave=1000,gdir="tmp",wait=TRUE) 
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
 if(mT1 > 0) modsT1 <- NULL
 if(mT1 == 0) modsT1 <- matrix(rep(0,s),nrow=1)
 for(i in 1:length(mc)) T1mods <- rbind(modsT1,modT1[[i]])
 
 j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
 t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
 t1size <- apply(T1mods,1,sum)
 T1mod <- data.frame(mod=t1mod,size=t1size)
 
 return(T1mod)
}

####

logBF0.fn <- function(msnps,t1pheno,G1.mat,t2pheno,G2.mat,trait1="T1",trait2="T2") {
 s <- length(msnps)
 nullmod <- matrix(c(rep(0,(2*s)),1),nrow=1) # no snp effects and only trait effect
 allsnps <- rbind(G1.mat[,msnps],G2.mat[t2pheno==1,msnps])
allpheno <- c(t1pheno,rep(2,sum(t2pheno)))
y <- character(length(allpheno))
y[allpheno==0] <- "CONTROL"
y[allpheno==1] <- trait1
y[allpheno==2] <- trait2
data1 <- data.frame(Y=y,allsnps)
m1 <- mlogit2logit(Y ~ 1|. -Y,data1,choices=c("CONTROL",trait1,trait2),base.choice=1)
Nullmod <- glib.1(x=m1$data[,(4+s):dim(m1$data)[2]],y=m1$data$Y.star,error="binomial", link="logit",models=nullmod,post.bymodel = FALSE)
logBF0 <- Nullmod$bf$twologB10[,1]*0.5
return(logBF0)
}

###

T1T2abfexp.fn <- function(T1bfexp,T2bfexp,logBF0) {
#' use to combine expanded model bfs from traits; all models

nT1 <- dim(T1bfexp)[1]
nT2 <- dim(T2bfexp)[1]
T1abf <- T1bfexp
T2abf <- T2bfexp

T1modsrep <- matrix(rep(t(T1abf),nT2),ncol=ncol(T1abf),byrow=TRUE)
T2modsrep <- matrix(rep(as.matrix(T2abf),each=nT1),ncol=ncol(T2abf),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

colnames(T1T2mods) <- c(paste("t1",names(T1abf),sep="."),paste("t2",names(T2abf),sep="."))
T1T2mods <- as.data.frame(T1T2mods,stringsAsFactors =FALSE)

T1T2mods$t1.lBF <- as.numeric(T1T2mods$t1.lBF)
T1T2mods$t2.lBF <- as.numeric(T1T2mods$t2.lBF)

t1t2 <- T1T2mods$t1.lBF+T1T2mods$t2.lBF+logBF0
out12 <- data.frame(t1t2.logbf=t1t2,T1T2mods)
names(out12) <- c("t1t2.logbf", "t1.str"  , "t1.tag"  ,   "t1.logbf" ,"t1.size"  ,  "t2.str" ,"t2.tag"   ,  "t2.logbf","t2.size")

out12$t1t2.Nexpmod <- 1
return(out12)
}

####

prior.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3) {
#'  called by PP.fn
  m1snp <- unlist(strsplit(as.character(bf$t1.str[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$t2.str[k]),"%"))
  m12snp <- intersect(m1snp,m2snp)
  nT1 <- as.numeric(bf$t1.size[k])
  nT2 <- as.numeric(bf$t2.size[k])
  nT1T2 <- length(m12snp)  
  prior12 <- shared*(nT1T2 > 0) + 1*(nT1T2 == 0)   # >1 if at least one overlap, 1 otherwise
  prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
  prior2 <- dbinom(nT2,size=s,prob=mT2/s)/choose(s,nT2)
  p <- prior1*prior2*prior12
  return(p)
  }

####

adjprior.fn <- function(prior,out) {
#'  called by PP.fn
 adjprior <- prior
 ind10 <- which(out$t2.size==0 & out$t1.size>0)
 ind01 <- which(out$t1.size==0 & out$t2.size>0)
 ind00 <- which(out$t2.size==0 & out$t1.size==0)
 ind11 <- which(out$t2.size>0 & out$t1.size>0)
 den <- sum(prior[ind11])
 den0 <- 0
 if(length(ind10) >0) den0 <- den0 + sum(prior[ind10])
 if(length(ind01) >0) den0 <- den0 + sum(prior[ind01])
 if(length(ind00) >0) den0 <- den0 + prior[ind00]
 adjprior[ind11] <- (1-den0)*prior[ind11]/den
 return(adjprior) 
}


####

PP.fn <- function(out,shared,s=length(snps(tags)),mT1=3,mT2=3,details=TRUE) {
  k.mat <- matrix(1:dim(out)[1],ncol=1)
  prior <- apply(k.mat,1,prior.fn,bf=out,shared=shared,s=s,mT1=mT1,mT2=mT2)
  if(min(out$t1.size)>0 & min(out$t2.size)>0 ) {
  adjprior <- prior/sum(prior*out$t1t2.Nexpmod)  
  } else { adjprior <- adjprior.fn(prior,out)
  }
  post <- exp(log(adjprior)+out[,"t1t2.logbf"]-max(out[,"t1t2.logbf"]))
  ppost <- post/sum(post*out$t1t2.Nexpmod)
  outPP <- matrix(ppost,nrow=1,dimnames=list(NULL,paste(out$t1.str,out$t2.str,sep="NEXT")))
 
 if(details) outPP <- data.frame(PP=ppost,prior=adjprior,logprior=log(adjprior),logbf=out[,"t1t2.logbf"],t1.mod=out$t1.str, t2.mod=out$t2.str,t1.logbf=out$t1.logbf,t2.logbf=out$t2.logbf,t1.size=out$t1.size,t2.size=out$t2.size)  
 
 return(outPP)
}


####

T1T2bfexp.low.fn <- function(BFkeep,tags,logBF0) {
#' at expanded marginal models with low PP (keep same bf as tag model) use counts of no. of expanded models to get bfs (and counts of these models) at  expanded joint models

tmp <- data.frame(logbf=BFkeep$t1.logbf,mod=BFkeep$t1.mod,size=BFkeep$t1.size)
nbest1 <- unique(tmp)
en1 <- expand.tags.bf(nbest1,tags)
T1lowexp <- en1[match(unique(en1$mod),en1$mod),c("logbf","mod","size")]
test<-table(en1$mod)
T1lowexp$Nexpmod <- test[match(T1lowexp$mod,names(test))]

tmp <- data.frame(logbf=BFkeep$t2.logbf,mod=BFkeep$t2.mod,size=BFkeep$t2.size)
nbest2 <- unique(tmp)
en2 <- expand.tags.bf(nbest2,tags)
T2lowexp <- en2[match(unique(en2$mod),en2$mod),c("logbf","mod","size")]
test<-table(en2$mod)
T2lowexp$Nexpmod <- test[match(T2lowexp$mod,names(test))]

b1 <- T1lowexp[match(BFkeep[,"t1.mod"],as.matrix(T1lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
b2 <- T2lowexp[match(BFkeep[,"t2.mod"],as.matrix(T2lowexp[,"mod"])),c("logbf","mod","size","Nexpmod")]
bf.all  <- data.frame(b1,b2)
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t1.Nexpmod","t2.logbf","t2.str", "t2.size","t2.Nexpmod")
    
t1t2bf <- bf.all$t1.logbf +bf.all$t2.logbf+logBF0
Nt1t2 <- bf.all$t1.Nexpmod *bf.all$t2.Nexpmod

bf.all$t1t2.logbf <- t1t2bf
bf.all$t1t2.Nexpmod <- Nt1t2

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
}

####

T1T2bfexp.best.fn <- function(T1T2abf,T1bfexp,T2bfexp,indbest,logBF0) {
#' on expanded marginal models, re-calculate joint bfs at "best" joint models so have joint bfs on expanded joint models

outkeep <- T1T2abf[indbest,]
names(T1bfexp) <- c("et1.logbf","t1.str","et1.size","et1.str","et1.rank")
b1 <- inner_join(outkeep,T1bfexp) # merge on t1.str to give all t1  expanded models 
names(T2bfexp) <- c("et2.logbf","t2.str","et2.size","et2.str","et2.rank")
b12 <- inner_join(b1,T2bfexp) # merge on t2.str

bf.all <- b12[,c("et1.logbf","et1.str", "et1.size","et2.logbf","et2.str", "et2.size")]
names(bf.all) <- c("t1.logbf","t1.str", "t1.size","t2.logbf","t2.str", "t2.size")

    
bf.all$t1t2.logbf <- bf.all$t1.logbf +bf.all$t2.logbf+logBF0

bf.all$t1t2.Nexpmod <- 1

out <- bf.all[order(bf.all$t1t2.logbf, decreasing = TRUE), ]
check0 <- which(out$t2.size==0)
if(length(check0)>0) out[check0,"t2.str"] <- 0
check0 <- which(out$t1.size==0)
if(length(check0)>0) out[check0,"t1.str"] <- 0

return(out)
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
#' @return a list with components pp (data.frame of joint PP, prior,
#'     logprior, logBF, and for each trait, model, logBF, model size)
#'     and logBF0, the offset term needed to calculate joint BFs
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

####
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


####

mod.split.fn <- function(k,PP,tnum) {
#' when two traits, tnum is 1 or 2 , i.e. 1st or 2nd trait
#' called by PP.marg.fn
tmp <- unlist(strsplit(row.names(PP)[k],"NEXT"))
msnp <- tmp[tnum]
if(is.na(msnp)) msnp <- "0" # null model for trait tnum
return(msnp)
}

####

PP.marg.fn <-function(PPall,K) {
 m<- dim(PPall)[2]
 PP <- PPall[,-m] # last column is no. of joint models with same PP 
PPm <- vector("list",K) 
 for(k in 1:K) {
  mod <- apply(matrix(1:dim(PP)[1],ncol=1),1,mod.split.fn,PP,k)
  tmp <- data.frame(mod,PP,Nmod=PPall[,m],row.names=NULL)
  pp <-  mergePP.fn(tmp)
  PPm[[k]] <- pp[order(pp[,1],decreasing=TRUE),]
    	}
          
return(PPm)
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
#' @title Calculate joint PP at a set of sharing scales
#' @param shared vector of sharing scales
#' @param bf output from BFexpanded.fn
#' @param tags set of tag SNPs for control set
#' @param trait1 name of trait 1, for filename purposes
#' @param trait2 name of trait 2, for filename purposes
#' @mydir directory to save output to
#' @return a list of PPmarg (PP for each trait) and mpp (marginal PP
#'     for each trait) evalauted at each sharing scale value
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
