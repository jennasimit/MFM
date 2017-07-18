
#msnps <- c(union(t1snps[,"var"],t2snps[,"var"]))
msnps <- ("rs2104286",  "rs56382813", "rs3118470" )

mT1 <- 0 # min model size for trait 1
mT2 <- 0 # min model size for trait 2
MT1 <- 3 # max model size for trait 1
MT2 <- 3 # max model size for trait 1
 
T1mods <- T1mods.fn(mT1,MT1,msnps)
T2mods <- T1mods.fn(mT2,MT2,msnps)

j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
t1size <- apply(T1mods,1,sum)
T1mod <- data.frame(mod=t1mod,size=t1size)
 
j.mat <- matrix(1:dim(T2mods)[1],ncol=1) 
t2mod <- apply(j.mat,1,format.mod.fn,T2mods)
t2size <- apply(T2mods,1,sum)
T2mod <- data.frame(mod=t2mod,size=t2size)
 
G1.mat <- as(G1,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)
 
G2.mat <- as(G2,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)
 

s<-244 # of tag snps in controls

t1size <- apply(T1mods,1,sum)
t2size <- apply(T2mods,1,sum)


bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 

bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 

bf1$size <- t1size
bf2$size <- t2size


prior.bin.fn <- function(k,bf,s,mT1=3) {
 #  s <- bf[k,"s"]
   #m1snp <- unlist(strsplit(as.character(bf$model[k]),"%")) 
    nT1 <- bf$size[k]
    prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
    return(prior1)
   }

 


M1 <- T1mods[-1,]
M2 <- T2mods[-1,]
ABF1 <- bf1[-1,"lBF"]
ABF2 <- bf2[-1,"lBF"]

p0 <- dbinom(0,s,3/s)

bf10 <- bf1[-1,]
bf20 <- bf2[-1,]
k.mat <- matrix(1:dim(bf10)[1],ncol=1)
pr1 <- apply(k.mat,1,prior.bin.fn,bf=bf10,s=s,mT1=3)
pr1 <- (1-p0)*pr1/sum(pr1)

pr2 <- apply(k.mat,1,prior.bin.fn,bf=bf20,s=s,mT1=3)
pr2 <- (1-p0)*pr1/sum(pr1)

pp12 <- marginalpp2(M1,M2,ABF1,ABF2,pr1,pr2,c(1,10,50,100),p0)


#################  previous way
mT1 <- 0
mT2 <- 0

 
T1mods <- T1mods.fn(mT1,MT1,msnps)
T2mods <- T1mods.fn(mT2,MT2,msnps)

j.mat <- matrix(1:dim(T1mods)[1],ncol=1) 
t1mod <- apply(j.mat,1,format.mod.fn,T1mods)
t1size <- apply(T1mods,1,sum)
T1mod <- data.frame(mod=t1mod,size=t1size)
 
j.mat <- matrix(1:dim(T2mods)[1],ncol=1) 
t2mod <- apply(j.mat,1,format.mod.fn,T2mods)
t2size <- apply(T2mods,1,sum)
T2mod <- data.frame(mod=t2mod,size=t2size)
 
G1.mat <- as(G1,"numeric")
data1 <- data.frame(Y=t1pheno,G1.mat)
 
G2.mat <- as(G2,"numeric")
data2 <- data.frame(Y=t2pheno,G2.mat)
 

s<-244

t1size <- apply(T1mods,1,sum)
t2size <- apply(T2mods,1,sum)


bf1 <- abf.calc(y=data1[,1],x=data1[,-1],models=T1mod$mod,family="binomial")[[1]] 

bf2 <- abf.calc(y=data2[,1],x=data2[,-1],models=T2mod$mod,family="binomial")[[1]] 

bf1$size <- t1size
bf2$size <- t2size


s <- length(msnps)
nullmod <- matrix(c(rep(0,(2*s)),1),nrow=1) # no snp effects and only trait effect
#t2m1snp <- t2snp.data[t2pheno==1,]
allsnps <- rbind(G1.mat[,msnps],G2.mat[t2pheno==1,msnps])
allpheno <- c(t1pheno,rep(2,sum(t2pheno)))
y <- character(length(allpheno))
y[allpheno==0] <- "CONTROL"
y[allpheno==1] <- trait1
y[allpheno==2] <- trait2
#G1.mat <- as(allsnps,"numeric")
data1 <- data.frame(Y=y,allsnps)
m1 <- mlogit2logit(Y ~ 1|. -Y,data1,choices=c("CONTROL",trait1,trait2),base.choice=1)
Nullmod <- glib.1(x=m1$data[,(4+s):dim(m1$data)[2]],y=m1$data$Y.star,error="binomial", link="logit",models=nullmod,post.bymodel = FALSE)
logBF0 <- Nullmod$bf$twologB10[,1]*0.5

T1T2abf <- T1T2abfexp.fn(bf1,bf2,logBF0)

shared <- c(1,10,50,100)
ns <- length(shared)
prior.fn <- function(k,bf,shared=1,s,mT1=3,mT2=3) {
#  s <- bf[k,"s"]
  m1snp <- unlist(strsplit(as.character(bf$t1.str[k]),"%"))
  m2snp <- unlist(strsplit(as.character(bf$t2.str[k]),"%"))
  m12snp <- intersect(m1snp,m2snp)
  nT1 <- as.numeric(bf$t1.size[k])
  nT2 <- as.numeric(bf$t2.size[k])
  nT1T2 <- length(m12snp)*(nT1!=0 | nT2!=0)
  #print(c(nT1,nT2,nT1T2))
  prior12 <- shared*(nT1T2 > 0) + 1*(nT1T2 == 0)   # >1 if at least one overlap, 1 otherwise
  #prior1 <- dbinom(nT1,size=s,prob=mT1/s)
  #prior2 <- dbinom(nT2,size=s,prob=mT2/s)
  prior1 <- dbinom(nT1,size=s,prob=mT1/s)/choose(s,nT1)
  prior2 <- dbinom(nT2,size=s,prob=mT2/s)/choose(s,nT2)
  p <- prior1*prior2*prior12
  return(p)
  }

adjprior.fn <- function(prior,out) {
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

##
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
 # indC <- order(colnames(tmp))
  #outPP <- tmp[,indC]
if(details) outPP <- data.frame(PP=ppost,prior=adjprior,logprior=log(adjprior),logbf=out[,"t1t2.logbf"],t1.mod=out$t1.str, t2.mod=out$t2.str,t1.logbf=out$t1.logbf,t2.logbf=out$t2.logbf,t1.size=out$t1.size,t2.size=out$t2.size)
  

  #ind <- order(ppost,decreasing=TRUE)
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1])  
  #outpp <- data.frame(ppost[ind],logABF=log10(exp(out[ind,1])),out[ind,-1]) 
  return(outPP)
}

PP <- NULL
for(k in 1:ns) PP <- rbind(PP,PP.fn(T1T2abf,shared=shared[k],s=244,mT1=3,mT2=3,details=FALSE))
row.names(PP) <- shared
PP <- t(PP)

PPall <- data.frame(PP,T1T2abf$t1t2.Nexpmod)


PP.marg.fn <- function(PPall,K) {
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

K=2
mpp <- vector("list",K)
PPmarg <- PP.marg.fn(PPall,K)


