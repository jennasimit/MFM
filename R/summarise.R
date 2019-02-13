#' @title internal function for MPP.fn
#' @param k index 
#' @param msep output from sep.fn
#' @param out storage vector for output
#' @param gnames SNP group names
#' @author Jenn Asimit
check.fn <- function(k,msep,out,gnames) {
     g <- length(gnames)
     p1 <- numeric(g) 
     for(j in 1:g) { 
     ind1 <- gnames[j] %in% msep[[k]]
     if(ind1) p1[j] <- out[k] 
     }
     return(p1)
    	}

#' @title internal function for MPP.fn
#' @param k index
#' @param mnames model names
#' @author Jenn Asimit
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


#' @title Marginal PP for models of a set of diseases, sharing information between the diseases
#' @param SM2 List of snpmod objects for a set of diseases
#' @param dis Vector of diseases for fine-mapping (subset from those in SM2)
#' @param thr Threshold such that the smallest set of models has cumulative PP >= thr
#' @param TOdds Vector of target odds of no sharing to sharing
#'	...
#' @param N0 number of shared controls
#' @param ND list of number of cases for a set of diseases
#' @param nsnps number of SNPs in the region
#' @return List consisting of PP: marginal PP for models and MPP: marginal PP of SNP inclusion
#' @export
#' @author Jenn Asimit
PPmarginal.multiple.fn <- function (SM2, dis, thr, TOdds, N0, ND,nsnps) 
{
	nd <- length(dis)    	
	kappas <- c()
	for(j in 1:length(TOdds)) kappas <- c(kappas,calckappa(nsnps=nsnps,p=2/nsnps,ndis=nd,target.odds=TOdds[j]))
    traits <- paste(dis, collapse = "-")
    bestmod.thr <- best.models(SM2[dis], cpp.thr = thr)
    M <- lapply(bestmod.thr, "[[", "str")
    pr <- lapply(bestmod.thr, "[[", "prior")
    abf <- lapply(bestmod.thr, "[[", "logABF")
    PP <- lapply(bestmod.thr, "[[", "PP")
    p0 <- snpprior(n = nsnps, expected = 2)["0"]
    STR = M[dis]
    ABF = abf[dis]
    PP <- PP[dis]
    pr = pr[dis]
    ND = ND[dis]
    message("\n\nCPP threshold = ", thr, "\n\tn.each (", paste(dis, 
        collapse = "/"), ") = ", paste(sapply(M[dis], length), 
        collapse = "/"))
    pp <- vector("list",length=nd) 
    #mpp <- vector("list",length=nd)  
       
     for(kappa in kappas) {
     ret <- marginalpp(STR, ABF, pr, kappa, p0, N0, ND,nsnps)    
     for(i in 1:nd) pp[[i]] <- cbind(pp[[i]],ret[[i]]$shared.pp)
     } 
      for(i in 1:nd) {
       pp[[i]] <- cbind(ret[[i]]$single.pp/sum(ret[[i]]$single.pp),pp[[i]])
       colnames(pp[[i]]) <- paste("pp",c("null",round(TOdds,2)),sep=".")
       rownames(pp[[i]]) <- rownames(ret[[i]])
     #  mpp[[i]] <- lapply(pp[[i]], MPP.fn)
       }
    mpp <- lapply(pp, MPP.fn)
    names(pp) <- dis
    mpp1 <- lapply(mpp, t)
    K <- length(dis)
    MPP <- mpp1[[1]] 
    for (k in 2:K) MPP <- smartbind(MPP, mpp1[[k]], fill = 0)
    return(list(PP = pp, MPP = MPP))
}


#' @title Approximate (fast) marginal PP for models of a set of diseases, sharing information between the diseases
#' @param SM2 List of snpmod objects for a set of diseases
#' @param dis Vector of diseases for fine-mapping (subset from those in SM2)
#' @param thr Threshold such that the smallest set of models has cumulative PP >= thr
#' @param kappa Vector of sharing values
#' @param fthr Second level filtering on all but first disease to give a faster approximation; retain SNPs within the smallest set of models such that cumulative PP >= thr
#' @param N0 number of shared controls
#' @param ND list of number of cases for a set of diseases
#' @param nsnps Number of SNPs in region
#' @return List consisting of PP: marginal PP for models and MPP: marginal PP of SNP inclusion
#' @export
#' @author Jenn Asimit
PPmarginal.multiple.fast.fn <- function(SM2,dis,thr,kappa,fthr,N0,ND,nsnps) {

traits <- paste(dis,collapse="-")
  
  bestmod.thr <- best.models(SM2[dis],cpp.thr=thr) 
  M <- lapply(bestmod.thr, "[[", "str") 
  pr <- lapply(bestmod.thr, "[[", "prior") 
  abf <- lapply(bestmod.thr, "[[", "logABF") 
  PP <- lapply(bestmod.thr, "[[", "PP") 
  p0 <- snpprior(n=nsnps,expected=2)["0"] 
  ND <- ND[dis]
  
  
 K <- length(dis) 

  
 pp <- vector("list",K)
 mpp <- vector("list",K)
 names(pp) <- dis
 names(mpp) <- dis
 for(j in 1:K) { # loop through each disease being first in vector to get its PP
  STR=M[dis] 
  ABF=abf[dis] 
  PP <- PP[dis] 
  pr=pr[dis]
  ND=ND[dis]
  message("\n\nCPP threshold = ",thr, "\n\tn.each (",paste(dis,collapse="/"),") = ",paste(sapply(M[dis],length),collapse="/")) 
  ret <- marginalone(STR,ABF,PP,pr,kappa,p0, fthr=fthr,N0,ND) # need to have each trait as 1st
  pp[[j]] <- ret$shared.pp
  rownames(pp[[j]]) <- ret$STR
  colnames(pp[[j]]) <- paste("pp",kappa,sep=".")
  mpp[[j]] <- MPP.fn(pp[[j]])
  dis <- c(dis[-1],dis[1])
  }

 MPP <- t(mpp[[1]])
for(k in 2:K) MPP <- smartbind(MPP, t(mpp[[k]]),fill=0)

return(list(PP=pp,MPP=MPP))
}


#' @title Filter probability matrix output by threshold 
#' @param pp Matrix of probabilities: rows are kappa values; columns are SNPs (MPP matrix) or models (PP matrix) 
#' @param pthr Threshold for keeping probabilities from SNPs/models that are greater than pthr for at least one kappa value
#' @return Probability matrix with the same rows as input matrix, but only columns with the maximum probability greater than the threshold are retained. 
pp.filter.fn <- function(pp,pthr) {
 maxp <- apply(pp,2,max)
 ind <- which(maxp > pthr)
 out <- pp[,ind]
 if(length(ind)==1) out <- matrix(pp[,ind],ncol=1,dimnames=list(rownames(pp),colnames(pp)[ind]))
 return(out)
 }


#' @title Plot a matrix of probabilities, where rows are kappa values and columns are SNPs (MPP matrix) or models (PP matrix)
#' @param pp Matrix of probabilities: rows are kappa values; columns are SNPs (MPP matrix) or models (PP matrix)
#' @param pthr Threshold for keeping probabilities from SNPs/models that are greater than pthr for at least one kappa value
#' @return Matrix of probabilites that are plotted (only SNPs/models with MPP/PP > pthr at some kappa value are retained)
#' @export
#' @author Jenn Asimit
PP.plot.fn <- function(pp,pthr) { 
 if(!is.na(pthr)) pp1 <- pp.filter.fn(pp,pthr=pthr)
 m <- dim(pp1)[1]
 p1 <- dim(pp1)[2]  
 
 ind <- grep("rs",colnames(pp1))
 tmp<-colnames(pp1)[ind]
 if(length(tmp)>0) {
  for(j in 1:length(tmp)) {
   tmp1 <- grep("rs",unlist(strsplit(tmp[j],"[.]")))
   colnames(pp1)[ind[j]] <- paste(unlist(strsplit(tmp[j],"[.]"))[-((tmp1+1):(tmp1+3))],collapse=".")
    }
	}
 	
 par(las=1) 
 #par(cex.axis=0.8) 
 #par(cex.axis=0.5) 
 image(1:m,1:p1,pp1,xlab="",ylab="",col=rev(heat.colors(128)),zlim=c(0,1.01),axes=FALSE)
 axis(2,at=1:p1,labels=colnames(pp1))
 box()
 axis(1,at=1:m,labels=rownames(pp1),las=3)
 image.plot(1:m,1:p1,pp1,xlab="",ylab="",col=rev(heat.colors(128)),zlim=c(0,1.01),legend.only=TRUE )
 return(pp1)
 		      	       		 }

#' @title SNP groups that have marginal PP of inclusion (MPPi) above a threshold for some disease
#' @param groups List consisting of groups: SNP groups formed for a set of diseases; summary: for each group the min r2, max r2 and MPPi for each disease is given; r2: matrix of r2 between each pair of SNPs
#' @param rdis Full set of diseases for fine-mapping 
#' @param mppi.thr Threshold for retaining SNP groups that have MPPi > mppi.thr for some disease in rdis
#' @return List of SNP groups that summarise fine-mapping results 
#' @export
make.snp.groups.fn <- function(groups,rdis,mppi.thr=.05) {
 group.summary <- groups$summary[,paste("sum.mppi.",rdis,sep="")]
 K <- length(rdis)
 check <- apply(group.summary,1,function(x) sum(x>=mppi.thr))
 keep <- which(check>0)
 snpGroups <- groups$groups@.Data[keep]
 names(snpGroups) <- LETTERS[1:length(snpGroups)] 
 return(snpGroups)
}


#' @title Summarise, for a set of diseases, posterior probabilities and marginal posterior probabilities by SNP groups
#' @param MPP Matrix of marginal probabilities for a set of diseases and kappa values; output as MPP from PPmarginal.multiple.fn
#' @param pp List consisting of posterior probability matrices for each disease; output as PP from PPmarginal.multiple.fn 
#' @param dis Vector of disease names
#' @param shared Vector of kappa values
#' @param snpGroups List of SNP groups for region
#' @return List consisting of mppGS: matrix of marginal probabilities by SNP groups; gPP: list of PP matrices by SNP groups
#' @export
#' @author Jenn Asimit
MPP.PP.groups.fn <- function(MPP,pp,dis,shared,snpGroups) {
 K <- length(dis)
alltraits <- dis

G <- character(length(colnames(MPP)))
ind <- which(colnames(MPP)=="1")
if(length(ind)>0) colnames(MPP)[ind] <- "m0"
for(k in 1:length(colnames(MPP))) {
check <- names(snpGroups)[grep(colnames(MPP)[k],snpGroups)];#print(c(k,colnames(MPP)[k],check))
if(length(check)!=0) {
G[k] <- check} else {G[k] <- colnames(MPP)[k]}
}

snpG <- unique(G)
ng <- length(snpGroups)
mppG <- matrix(0,ncol=ng,nrow=dim(MPP)[1],dimnames=list(rownames(MPP),names(snpGroups)))
for(k in 1:ng) mppG[,k] <- apply(as.matrix(MPP[,which(G==names(snpGroups)[k])]),1,sum)

mppGS <- mppG
notinG <- setdiff(snpG,names(snpGroups))
notinG <- setdiff(notinG,"m0")
if(length(notinG)>0) {
mppS <- matrix(0,ncol=length(notinG),nrow=dim(MPP)[1],dimnames=list(rownames(MPP),notinG))
for(k in 1:length(notinG)) mppS[,k] <- apply(as.matrix(MPP[,which(G==notinG[k])]),1,sum)
mppGS <- cbind(mppG,mppS)
}
rownames(mppGS) <- paste(rep(alltraits,each=length(shared)),shared,sep=".")


gPP <- vector("list",K)
cmpp <- colnames(MPP)

# mppGS output
# find pp and pp by group

PPmarg <- pp
for(k in 1:K) {
ind <- which(rownames(PPmarg[[k]])=="1")
if (length(ind) > 0) {
            rownames(PPmarg[[k]])[ind] <- "null"
        }
	tmp <- rownames(PPmarg[[k]])
        rn <- character(length(tmp))
        Grn <- character(length(tmp))
        for (l in 1:length(tmp)) {
            tmp1 <- unlist(strsplit(tmp[l], "%"))
			rn[l] <- tmp[l]
            gg <- NULL
            for (ll in 1:length(tmp1)) gg <- c(gg, G[grep(tmp1[ll],colnames(MPP))])
            Grn[l] <- paste(gg[order(gg)], collapse = ".")	
            if(Grn[l] == "") Grn[l] <- "null"
        }
        rownames(PPmarg[[k]]) <- rn
        gPPmarg <- PPmarg[[k]]
        mods <- unique(Grn)
        nm <- length(mods)
        gPP[[k]] <- matrix(0, ncol = length(shared), nrow = nm,
            dimnames = list(mods, shared))
        for (mm in 1:nm) {
         mat <- as.matrix(gPPmarg[Grn == mods[mm], ], ncol = length(shared), byrow = FALSE)
         if(dim(mat)[2]==1) {
                gPP[[k]][mm, ] <- as.vector(mat)
                } else { gPP[[k]][mm, ] <- apply(mat, 2, sum) }
         }
    }
    return(list(mppGS = mppGS, gPP = gPP))

}


