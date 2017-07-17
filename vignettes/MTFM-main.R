
# load traits 1 and 2
# t1snp.data, t2snp.data: SnpMatrix objects for traits 1 and 2 
# t1pheno, t2pheno: phenotype vectors of 0,1 for controls,cases of traits 1 and 2; have common set of controls here

# trait1, trait2 are names of traits
trait1="T1"
trait2="T2"

# specify r2 threshold for defining tag SNPs 
r2=0.9

# specify "best" SNPs threshold such that tag SNPs with marginal PP > mmppthr are carried forward in models
mppthr=0.001

# specify minimum and maximum model sizes for traits 1 and 2
mT1=0; MT1=3
mT2=0; MT2=3

# specify directory to save files
mydir <- "tmp"

# specify sharing scales
shared <- c(1,5,10,20,100,1000)

tags <- make.tags.fn(r2=r2,G0=t1snp.data[t1pheno==0,],mydir=NA) # set mydir=NA if do not want to save tags to an Rdata file

# find joint BFs and then posterior probabilities for tag SNP models
ppbf1 <- PPBF.tags.fn(t1snp.data,t1pheno,t2snp.data,t2pheno,mydir,tags,mppthr,mT1,MT1,mT2,MT2,trait1,trait2)

# expand tag SNP models 
bf <- BFexpanded.fn(ppbf1,tags,t1snp.data,t1pheno,t2snp.data,t2pheno,trait1,trait2)


PPfinal <- PPshared.fn(shared,bf,tags,trait1,trait2,mydir)